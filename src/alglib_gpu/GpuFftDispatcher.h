#pragma once

#include "GpuBackend.h"

#include <atomic>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

namespace alglib_gpu
{
class GpuFftDispatcher
  {
public:
   GpuFftDispatcher();
   ~GpuFftDispatcher();

   JobStatus initialize(const GpuConfig& cfg);
   void      shutdown();

   JobStatus submit(const double* data,
                    int           fft_len,
                    FftKind       kind,
                    bool          inverse,
                    std::uint64_t& handle_out);

   JobStatus collect(std::uint64_t handle, double* buffer, int buffer_len, bool& ready);

   JobStatus status(std::uint64_t handle, bool& ready) const;

private:
   struct FftJob
     {
      std::uint64_t       handle;
      FftKind             kind;
      bool                inverse;
      int                 fft_len;
      std::vector<double> input;
      std::vector<double> output;
      JobStatus           status;
      bool                completed;
     };

   void worker_loop();

private:
   BackendContext                                              backend_;
   GpuConfig                                                   config_;
   std::vector<std::thread>                                    workers_;
   std::unordered_map<std::uint64_t, std::shared_ptr<FftJob>>  jobs_;
   std::deque<std::shared_ptr<FftJob>>                         queue_;
   mutable std::mutex                                          jobs_mutex_;
   std::mutex                                                  queue_mutex_;
   std::condition_variable                                     queue_cv_;
   std::atomic<uint64_t>                                       next_handle_;
   std::atomic<bool>                                           running_;
  };
} // namespace alglib_gpu
