#include "GpuFftDispatcher.h"

#include <algorithm>

namespace alglib_gpu
{
GpuFftDispatcher::GpuFftDispatcher()
    : next_handle_(1),
      running_(false)
  {
  }

GpuFftDispatcher::~GpuFftDispatcher()
  {
   shutdown();
  }

JobStatus GpuFftDispatcher::initialize(const GpuConfig& cfg)
  {
   if(running_)
      return JobStatus::Ok;

   config_ = cfg;
   auto status = backend_.initialize(cfg);
   if(status != JobStatus::Ok)
      return status;

   running_ = true;

   const int worker_count = std::max(1, cfg.stream_count);
   workers_.reserve(worker_count);
   for(int i = 0; i < worker_count; ++i)
      workers_.emplace_back([this]() { worker_loop(); });

   return JobStatus::Ok;
  }

void GpuFftDispatcher::shutdown()
  {
   running_ = false;
   queue_cv_.notify_all();

   for(auto& worker : workers_)
     {
      if(worker.joinable())
         worker.join();
     }
   workers_.clear();

   {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    queue_.clear();
   }

   {
    std::lock_guard<std::mutex> lock(jobs_mutex_);
    jobs_.clear();
   }

   backend_.shutdown();
  }

JobStatus GpuFftDispatcher::submit(const double* data,
                                   int           fft_len,
                                   FftKind       kind,
                                   bool          inverse,
                                   std::uint64_t& handle_out)
  {
   if(data == nullptr || fft_len <= 0)
      return JobStatus::ErrorInvalidArgs;

   auto job          = std::make_shared<FftJob>();
   job->handle       = next_handle_.fetch_add(1);
   job->kind         = kind;
   job->inverse      = inverse;
   job->fft_len      = fft_len;
   job->status       = JobStatus::Pending;
   job->completed    = false;

   const int expected_len = (kind == FftKind::Complex)
                                ? fft_len * 2
                                : (inverse ? fft_len * 2 : fft_len);
   job->input.assign(data, data + expected_len);

   {
    std::lock_guard<std::mutex> lock(jobs_mutex_);
    jobs_.emplace(job->handle, job);
   }
   {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    queue_.push_back(std::move(job));
   }
   queue_cv_.notify_one();

   handle_out = job->handle;
   return JobStatus::Pending;
  }

JobStatus GpuFftDispatcher::collect(std::uint64_t handle, double* buffer, int buffer_len, bool& ready)
  {
   std::shared_ptr<FftJob> job;
   {
    std::lock_guard<std::mutex> lock(jobs_mutex_);
    auto                        it = jobs_.find(handle);
    if(it == jobs_.end())
       return JobStatus::ErrorInvalidArgs;
    job = it->second;
   }

   ready = job->completed;
   if(!job->completed)
      return JobStatus::Pending;

   const int expected = static_cast<int>(job->output.size());
   if(buffer == nullptr || buffer_len < expected)
      return JobStatus::ErrorInvalidArgs;

   std::copy(job->output.begin(), job->output.end(), buffer);

   std::lock_guard<std::mutex> lock(jobs_mutex_);
   jobs_.erase(handle);

   return job->status;
  }

JobStatus GpuFftDispatcher::status(std::uint64_t handle, bool& ready) const
  {
   std::lock_guard<std::mutex> lock(jobs_mutex_);
   auto                        it = jobs_.find(handle);
   if(it == jobs_.end())
      return JobStatus::ErrorInvalidArgs;
   ready = it->second->completed;
   return it->second->status;
  }

void GpuFftDispatcher::worker_loop()
  {
   while(running_)
     {
      std::shared_ptr<FftJob> job;
      {
       std::unique_lock<std::mutex> lock(queue_mutex_);
       queue_cv_.wait(lock, [this]() { return !queue_.empty() || !running_; });
       if(!running_ && queue_.empty())
          return;

       job = queue_.front();
       queue_.pop_front();
      }

      if(!job)
         continue;

      FftSubmitDesc desc{};
      desc.data    = job->input.data();
      desc.fft_len = job->fft_len;
      desc.kind    = job->kind;
      desc.inverse = job->inverse;

      std::vector<double> buffer;
      const auto          status = backend_.execute_fft(desc, buffer);

      job->status    = status;
      job->completed = true;
      job->output    = std::move(buffer);
     }
  }
} // namespace alglib_gpu
