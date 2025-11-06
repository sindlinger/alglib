#pragma once

#include "alglib_pipe_messages.h"
#include "GpuBackend.h"
#include "dataanalysis.h"

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gpu_service
{
// Legacy FFT payload removido – usar apenas PROCESS_* para FFT

struct ProcessResultPayload
  {
   wave_pipe::Status status;
   std::vector<double> primary;
   std::vector<double> secondary;
   std::vector<wave_pipe::FftCyclePayload> cycles;
   std::vector<wave_pipe::InstantMetricsPayload> instant;
   std::vector<wave_pipe::TransitionPayload> transitions;
   std::vector<std::uint8_t> metadata;
  };

struct BatchResultPayload
  {
   wave_pipe::Status status;
   std::vector<std::uint8_t> buffer;
  };

class ServiceRuntime
  {
public:
   struct Metadata
     {
      std::string version;
      std::string backend;
      std::string driver;
      std::string build_timestamp;
      std::string host;
      int         workers = 0;
      std::uint64_t pid = 0;
     };

   ServiceRuntime();
   ~ServiceRuntime();

   bool start(const alglib_gpu::GpuConfig& cfg, int worker_count);
   void stop();

   void set_metadata(const Metadata& meta);

   // FFT via PROCESS_* (sem submit_fft/fetch_fft legados)

   void                submit_process(const wave_pipe::ProcessRequest& header,
                                      const double* primary,
                                      const double* secondary,
                                      const std::uint8_t* params);
   wave_pipe::Status   fetch_process(std::uint64_t tag, ProcessResultPayload& out);

   void                submit_batch(const wave_pipe::BatchHeader& header,
                                    const wave_pipe::BatchTaskDescriptor* tasks,
                                    const std::uint8_t* data_blob,
                                    const std::uint8_t* params_blob);
   wave_pipe::Status   fetch_batch(std::uint64_t tag, BatchResultPayload& out);

   alglib_gpu::BackendType resolved_backend() const;

private:
   enum class JobType { Process, Batch };

   // FftJob removido

   struct ProcessJob
     {
      wave_pipe::ProcessRequest header;
      std::vector<double>       primary;
      std::vector<double>       secondary;
      std::vector<std::uint8_t> params;
     };

  struct BatchJob
    {
     wave_pipe::BatchHeader header;
     std::vector<wave_pipe::BatchTaskDescriptor> tasks;
     std::vector<std::uint8_t> data_blob;
     std::vector<std::uint8_t> params_blob;
    };

   struct SsaContext
    {
     std::unique_ptr<alglib::ssamodel> model;
     int window          = 0;
     int topk_direct     = 0;
     int topk_realtime   = 0;
    };

   struct JobWrapper
     {
      JobType type;
      std::shared_ptr<ProcessJob> process;
      std::shared_ptr<BatchJob>   batch;
     };

   void worker_loop();
   void process_signal(const std::shared_ptr<ProcessJob>& job);
   void process_batch(const std::shared_ptr<BatchJob>& job);

   static wave_pipe::Status convert_status(alglib_gpu::JobStatus status);

private:
   alglib_gpu::BackendContext       backend_;
   alglib_gpu::GpuConfig            config_;

   std::vector<std::thread>         workers_;
   std::deque<JobWrapper>           queue_;
   std::mutex                       queue_mutex_;
   std::condition_variable          queue_cv_;
   bool                             running_;

   std::mutex                       state_mutex_;
   // pendências FFT legadas removidas
   std::unordered_set<std::uint64_t> pending_process_;
   std::unordered_map<std::uint64_t, ProcessResultPayload> process_results_;
   std::unordered_set<std::uint64_t> pending_batch_;
   std::unordered_map<std::uint64_t, BatchResultPayload>   batch_results_;

   Metadata                      metadata_;

   std::mutex                       ssa_mutex_;
   std::unordered_map<std::uint64_t, SsaContext> ssa_contexts_;
   std::atomic<std::uint64_t>       next_ssa_handle_{1};
  };
} // namespace gpu_service
