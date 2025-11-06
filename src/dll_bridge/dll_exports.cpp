#include "ap.h"
#include "dataanalysis.h"
#include "fasttransforms.h"
#include "GpuCommon.h"
#include "GpuFftDispatcher.h"
#include "alglib_pipe_messages.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#undef ERROR
#undef OK
#undef INVALID
#undef PENDING
#endif
#include <mutex>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace alglib;

namespace
{
std::mutex                         g_fft_mutex;
std::atomic<bool>                  g_fft_ready{false};
alglib_gpu::GpuFftDispatcher       g_fft_dispatcher; // legacy local FFT path (deprecated)

#ifdef _WIN32
HANDLE                             g_pipe_handle = INVALID_HANDLE_VALUE;
std::mutex                         g_pipe_mutex;
std::atomic<std::uint64_t>         g_next_pipe_tag{1};
#endif

constexpr int kDefaultStreamCount          = 8;
constexpr int kDefaultMinWindow            = 1024;
constexpr int kDefaultJobTimeoutMs         = 5000;

std::atomic<int> g_active_timeout_ms{kDefaultJobTimeoutMs};

constexpr int kCycleStride      = 5;
constexpr int kInstantStride    = 6;
constexpr int kTransitionStride = 5;

struct ProcessResult
  {
   wave_pipe::Status status{wave_pipe::Status::ERROR};
   std::vector<double> primary;
   std::vector<double> secondary;
   std::vector<wave_pipe::FftCyclePayload> cycles;
   std::vector<wave_pipe::InstantMetricsPayload> instant;
   std::vector<wave_pipe::TransitionPayload> transitions;
  };

struct BatchProcessBuffer
  {
   wave_pipe::Status status{wave_pipe::Status::ERROR};
   std::vector<std::uint8_t> payload;
  };

std::mutex                                       g_process_async_mutex;
std::unordered_set<std::uint64_t>                g_process_pending_handles;
std::unordered_map<std::uint64_t, ProcessResult> g_process_ready_results;

std::mutex                                       g_ssa_proxy_mutex;
std::unordered_map<alglib::ssamodel*, std::uint64_t> g_ssa_proxy_handles;

struct ResampleParams
  {
   double factor;
   double cutoff;
   std::int32_t method;
   std::int32_t reserved;
  };

struct ZeroPadParams
  {
   std::int32_t pad_left;
   std::int32_t pad_right;
   std::int32_t reserved1;
   std::int32_t reserved2;
  };

struct RemoveDcParams
  {
   double alpha;
   std::int32_t mode;
   std::int32_t reserved1;
   std::int32_t reserved2;
  };

struct GaussianParams
  {
   double sigma_low;
   double sigma_high;
   double gain;
   std::int32_t mode;
  };

struct NotchParams
  {
   double center;
   double width;
   double depth;
   double gain;
  };

struct SpectrumParams
  {
   std::int32_t top_k;
   double       min_period;
   double       max_period;
   double       energy_threshold;
   double       upscale_factor;
  };

struct TransitionParams
  {
   double       energy_threshold;
   double       phase_jump_threshold;
   std::int32_t min_duration;
   std::int32_t type_mask;
  };

struct InstantParams
  {
   std::int32_t smooth_window;
   double       epsilon;
  };

struct MaskParams
  {
   std::int32_t mode;
   std::int32_t reserved1;
   std::int32_t reserved2;
   std::int32_t reserved3;
  };

struct DenoiseParams
  {
   std::int32_t method;
   double       threshold;
   double       beta;
   std::int32_t iterations;
  };

struct UpscaleParams
  {
   double       factor;
   std::int32_t mode;
   std::int32_t normalize;
   std::int32_t reserved;
  };

struct DownscaleParams
  {
   double       factor;
   std::int32_t mode;
   std::int32_t anti_alias;
   std::int32_t reserved;
  };

struct ConvolutionParams
  {
   std::int32_t normalize;
   std::int32_t reserved1;
   std::int32_t reserved2;
   std::int32_t reserved3;
  };

struct PhaseParams
  {
   std::int32_t method;
   std::int32_t reserved1;
   std::int32_t reserved2;
   std::int32_t reserved3;
  };

struct SsaCreateParams
  {
   std::int32_t window;
   std::int32_t topk_direct;
   std::int32_t topk_realtime;
   std::int32_t reserved;
  };

struct SsaHandleParam
  {
   double handle;
   double value;
   double extra1;
   double extra2;
  };

struct SsaConfigParams
  {
   double       handle;
   std::int32_t window;
   std::int32_t topk_direct;
   std::int32_t topk_realtime;
  };

struct SsaAppendParams
  {
   double handle;
   double value;
   double update_iters;
   double reserved;
  };

struct SsaAnalyzeParams
  {
   double       handle;
   std::int32_t nticks;
   std::int32_t reserved1;
   std::int32_t reserved2;
  };

struct SsaForecastParams
  {
   double       handle;
   std::int32_t nticks;
   std::int32_t reserved1;
   std::int32_t reserved2;
  };

void store_ready_result(std::uint64_t tag, ProcessResult&& result)
  {
   std::lock_guard<std::mutex> lock(g_process_async_mutex);
   g_process_ready_results[tag] = std::move(result);
   g_process_pending_handles.erase(tag);
  }

void mark_handle_pending(std::uint64_t tag)
  {
   std::lock_guard<std::mutex> lock(g_process_async_mutex);
   g_process_pending_handles.insert(tag);
  }

bool take_ready_result(std::uint64_t tag, ProcessResult& out)
  {
   std::lock_guard<std::mutex> lock(g_process_async_mutex);
   auto it = g_process_ready_results.find(tag);
   if(it == g_process_ready_results.end())
      return false;
   out = std::move(it->second);
   g_process_ready_results.erase(it);
   g_process_pending_handles.erase(tag);
   return true;
  }

bool handle_is_pending(std::uint64_t tag)
  {
   std::lock_guard<std::mutex> lock(g_process_async_mutex);
   return g_process_pending_handles.find(tag) != g_process_pending_handles.end();
  }

void clear_pending_handle(std::uint64_t tag)
  {
   std::lock_guard<std::mutex> lock(g_process_async_mutex);
   g_process_pending_handles.erase(tag);
   g_process_ready_results.erase(tag);
  }

alglib_gpu::JobStatus copy_primary_output(const ProcessResult& result,
                                          double* output,
                                          int     capacity,
                                          int*    out_len)
  {
   const int produced = static_cast<int>(result.primary.size());
   if(out_len)
      *out_len = produced;

   if(produced == 0)
      return alglib_gpu::JobStatus::Ok;

   if(output == nullptr || capacity < produced)
      return alglib_gpu::JobStatus::ErrorInvalidArgs;

   std::copy(result.primary.begin(), result.primary.end(), output);
  return alglib_gpu::JobStatus::Ok;
 }

std::uint64_t decode_handle(double value)
  {
   if(!std::isfinite(value) || value < 0.0)
      return 0;
   return static_cast<std::uint64_t>(std::llround(value));
  }

bool resolve_ssa_handle(alglib::ssamodel* model, std::uint64_t& handle)
  {
   std::lock_guard<std::mutex> lock(g_ssa_proxy_mutex);
   auto it = g_ssa_proxy_handles.find(model);
   if(it == g_ssa_proxy_handles.end())
      return false;
   handle = it->second;
   return true;
  }

void register_ssa_handle(alglib::ssamodel* model, std::uint64_t handle)
  {
   std::lock_guard<std::mutex> lock(g_ssa_proxy_mutex);
   g_ssa_proxy_handles[model] = handle;
  }

void unregister_ssa_handle(alglib::ssamodel* model)
  {
   std::lock_guard<std::mutex> lock(g_ssa_proxy_mutex);
   g_ssa_proxy_handles.erase(model);
  }

alglib_gpu::JobStatus copy_secondary_output(const ProcessResult& result,
                                            double* output,
                                            int     capacity,
                                            int*    out_len)
  {
   const int produced = static_cast<int>(result.secondary.size());
   if(out_len)
      *out_len = produced;

   if(produced == 0)
      return alglib_gpu::JobStatus::Ok;

   if(output == nullptr || capacity < produced)
      return alglib_gpu::JobStatus::ErrorInvalidArgs;

   std::copy(result.secondary.begin(), result.secondary.end(), output);
   return alglib_gpu::JobStatus::Ok;
  }

alglib_gpu::JobStatus copy_cycles_output(const ProcessResult& result,
                                         double* cycles_out,
                                         int     stride,
                                         int     capacity,
                                         int*    cycles_len)
  {
   const int count = static_cast<int>(result.cycles.size());
   if(cycles_len)
      *cycles_len = count;

   if(count == 0)
      return alglib_gpu::JobStatus::Ok;

   if(cycles_out == nullptr || stride < kCycleStride || capacity < count)
      return alglib_gpu::JobStatus::ErrorInvalidArgs;

   for(int i = 0; i < count; ++i)
     {
      const auto& cycle = result.cycles[static_cast<std::size_t>(i)];
      double* row = cycles_out + static_cast<std::size_t>(i) * stride;
      row[0] = static_cast<double>(cycle.index);
      row[1] = cycle.amplitude;
      row[2] = cycle.phase;
      row[3] = cycle.frequency;
      row[4] = cycle.period;
     }
   return alglib_gpu::JobStatus::Ok;
  }

alglib_gpu::JobStatus copy_instant_output(const ProcessResult& result,
                                          double* instant_out,
                                          int     stride,
                                          int     capacity,
                                          int*    instant_len)
  {
   const int count = static_cast<int>(result.instant.size());
   if(instant_len)
      *instant_len = count;

   if(count == 0)
      return alglib_gpu::JobStatus::Ok;

   if(instant_out == nullptr || stride < kInstantStride || capacity < count)
      return alglib_gpu::JobStatus::ErrorInvalidArgs;

   for(int i = 0; i < count; ++i)
     {
      const auto& inst = result.instant[static_cast<std::size_t>(i)];
      double* row = instant_out + static_cast<std::size_t>(i) * stride;
      row[0] = inst.amplitude;
      row[1] = inst.frequency;
      row[2] = inst.energy;
      row[3] = inst.phase;
      row[4] = inst.envelope_up;
      row[5] = inst.envelope_down;
     }
   return alglib_gpu::JobStatus::Ok;
  }

alglib_gpu::JobStatus copy_transitions_output(const ProcessResult& result,
                                              double* transitions_out,
                                              int     stride,
                                              int     capacity,
                                              int*    transition_len)
  {
   const int count = static_cast<int>(result.transitions.size());
   if(transition_len)
      *transition_len = count;

   if(count == 0)
      return alglib_gpu::JobStatus::Ok;

   if(transitions_out == nullptr || stride < kTransitionStride || capacity < count)
      return alglib_gpu::JobStatus::ErrorInvalidArgs;

   for(int i = 0; i < count; ++i)
     {
      const auto& tr = result.transitions[static_cast<std::size_t>(i)];
      double* row = transitions_out + static_cast<std::size_t>(i) * stride;
      row[0] = static_cast<double>(tr.start_index);
      row[1] = static_cast<double>(tr.end_index);
      row[2] = tr.energy_ratio;
      row[3] = tr.phase_shift;
      row[4] = static_cast<double>(tr.type);
     }
   return alglib_gpu::JobStatus::Ok;
  }

#ifdef _WIN32
constexpr wchar_t kDefaultPipeName[] = L"\\\\.\\pipe\\alglib-wave_pipe";

HANDLE ensure_pipe_connected()
  {
   if(g_pipe_handle != INVALID_HANDLE_VALUE)
      return g_pipe_handle;

   for(int attempt = 0; attempt < 3; ++attempt)
     {
      if(WaitNamedPipeW(kDefaultPipeName, 1000) || GetLastError() == ERROR_SUCCESS)
        {
         HANDLE handle = CreateFileW(kDefaultPipeName,
                                     GENERIC_READ | GENERIC_WRITE,
                                     0,
                                     nullptr,
                                     OPEN_EXISTING,
                                     0,
                                     nullptr);
         if(handle != INVALID_HANDLE_VALUE)
           {
            DWORD mode = PIPE_READMODE_MESSAGE;
            SetNamedPipeHandleState(handle, &mode, nullptr, nullptr);
            g_pipe_handle = handle;
            return g_pipe_handle;
           }
        }
      Sleep(100);
     }

   return INVALID_HANDLE_VALUE;
  }

void close_pipe()
  {
   if(g_pipe_handle != INVALID_HANDLE_VALUE)
     {
      CloseHandle(g_pipe_handle);
      g_pipe_handle = INVALID_HANDLE_VALUE;
     }
  }

bool pipe_write(HANDLE pipe, const void* data, size_t bytes)
  {
   const std::uint8_t* ptr = reinterpret_cast<const std::uint8_t*>(data);
   size_t remaining = bytes;
   while(remaining > 0)
     {
      DWORD written = 0;
      if(!WriteFile(pipe, ptr, static_cast<DWORD>(remaining), &written, nullptr) || written == 0)
         return false;
      remaining -= written;
      ptr += written;
     }
   return true;
  }

bool pipe_read(HANDLE pipe, void* data, size_t bytes)
  {
   std::uint8_t* ptr = reinterpret_cast<std::uint8_t*>(data);
   size_t        remaining = bytes;
   while(remaining > 0)
     {
      DWORD read = 0;
      if(!ReadFile(pipe, ptr, static_cast<DWORD>(remaining), &read, nullptr) || read == 0)
         return false;
      remaining -= read;
      ptr += read;
     }
   return true;
  }

alglib_gpu::JobStatus status_from_pipe(wave_pipe::Status status)
  {
   switch(status)
     {
      case wave_pipe::Status::OK:
         return alglib_gpu::JobStatus::Ok;
      case wave_pipe::Status::PENDING:
         return alglib_gpu::JobStatus::Pending;
      case wave_pipe::Status::TIMEOUT:
         return alglib_gpu::JobStatus::ErrorTimeout;
      case wave_pipe::Status::INVALID:
         return alglib_gpu::JobStatus::ErrorInvalidArgs;
      case wave_pipe::Status::ERROR:
      default:
         return alglib_gpu::JobStatus::ErrorBackendFailure;
     }
  }

alglib_gpu::JobStatus submit_process_async(wave_pipe::Operation operation,
                                           const double* primary,
                                           int           primary_len,
                                           const double* secondary,
                                           int           secondary_len,
                                           const void*   params,
                                           std::uint32_t param_size,
                                           std::uint64_t& handle_out)
  {
   handle_out = 0;

   if(primary == nullptr || primary_len <= 0)
      return alglib_gpu::JobStatus::ErrorInvalidArgs;

   std::lock_guard<std::mutex> lock(g_pipe_mutex);

   HANDLE pipe = ensure_pipe_connected();
   if(pipe == INVALID_HANDLE_VALUE)
      return alglib_gpu::JobStatus::ErrorNoBackend;

   const std::uint64_t tag = g_next_pipe_tag.fetch_add(1);

   wave_pipe::ProcessRequest req{};
   req.magic      = wave_pipe::MESSAGE_MAGIC;
   req.version    = wave_pipe::PROTOCOL_VERSION;
   req.cmd        = static_cast<std::uint32_t>(wave_pipe::Command::PROCESS_SUBMIT);
   req.user_tag   = tag;
   req.operation  = static_cast<std::uint32_t>(operation);
   req.flags      = 0;
   req.window_len = primary_len;
   req.aux_len    = secondary_len;
   req.param_size = param_size;

   if(!pipe_write(pipe, &req, sizeof(req)))
     {
      close_pipe();
      return alglib_gpu::JobStatus::ErrorBackendFailure;
     }

   const size_t primary_bytes = static_cast<size_t>(primary_len) * sizeof(double);
   if(primary_bytes > 0)
     {
      if(!pipe_write(pipe, primary, primary_bytes))
        {
         close_pipe();
         return alglib_gpu::JobStatus::ErrorBackendFailure;
        }
     }

   const size_t secondary_bytes = (secondary && secondary_len > 0)
                                      ? static_cast<size_t>(secondary_len) * sizeof(double)
                                      : 0;
   if(secondary_bytes > 0)
     {
      if(!pipe_write(pipe, secondary, secondary_bytes))
        {
         close_pipe();
         return alglib_gpu::JobStatus::ErrorBackendFailure;
        }
     }

   if(param_size > 0 && params != nullptr)
     {
      if(!pipe_write(pipe, params, param_size))
        {
         close_pipe();
         return alglib_gpu::JobStatus::ErrorBackendFailure;
        }
     }

   FlushFileBuffers(pipe);

   wave_pipe::ProcessResponse response{};
   if(!pipe_read(pipe, &response, sizeof(response)))
     {
      close_pipe();
      return alglib_gpu::JobStatus::ErrorBackendFailure;
     }

   ProcessResult local{};
   local.status = static_cast<wave_pipe::Status>(response.status);

   if(local.status == wave_pipe::Status::OK)
     {
      if(response.primary_count > 0)
        {
         local.primary.resize(response.primary_count);
         const size_t bytes = static_cast<size_t>(response.primary_count) * sizeof(double);
         if(!pipe_read(pipe, local.primary.data(), bytes))
           {
            close_pipe();
            return alglib_gpu::JobStatus::ErrorBackendFailure;
           }
        }

      if(response.secondary_count > 0)
        {
         local.secondary.resize(response.secondary_count);
         const size_t bytes = static_cast<size_t>(response.secondary_count) * sizeof(double);
         if(!pipe_read(pipe, local.secondary.data(), bytes))
           {
            close_pipe();
            return alglib_gpu::JobStatus::ErrorBackendFailure;
           }
        }

      if(response.cycle_count > 0)
        {
         local.cycles.resize(response.cycle_count);
         const size_t bytes = static_cast<size_t>(response.cycle_count) * sizeof(wave_pipe::FftCyclePayload);
         if(!pipe_read(pipe, local.cycles.data(), bytes))
           {
            close_pipe();
            return alglib_gpu::JobStatus::ErrorBackendFailure;
           }
        }

      if(response.instant_count > 0)
        {
         local.instant.resize(response.instant_count);
         const size_t bytes = static_cast<size_t>(response.instant_count) * sizeof(wave_pipe::InstantMetricsPayload);
         if(!pipe_read(pipe, local.instant.data(), bytes))
           {
            close_pipe();
            return alglib_gpu::JobStatus::ErrorBackendFailure;
           }
        }

      if(response.transition_count > 0)
        {
         local.transitions.resize(response.transition_count);
         const size_t bytes = static_cast<size_t>(response.transition_count) * sizeof(wave_pipe::TransitionPayload);
         if(!pipe_read(pipe, local.transitions.data(), bytes))
           {
            close_pipe();
            return alglib_gpu::JobStatus::ErrorBackendFailure;
           }
        }
     }

   handle_out = tag;

   if(local.status == wave_pipe::Status::PENDING)
     {
      mark_handle_pending(tag);
     }
   else
     {
      store_ready_result(tag, std::move(local));
     }

   return status_from_pipe(local.status);
  }

alglib_gpu::JobStatus fetch_process_async(std::uint64_t tag,
                                          ProcessResult& result,
                                          int           wait_ms)
  {
   result = ProcessResult{};

   if(take_ready_result(tag, result))
      return status_from_pipe(result.status);

   if(!handle_is_pending(tag))
      return alglib_gpu::JobStatus::ErrorInvalidArgs;

   const bool has_timeout = wait_ms > 0;
   const auto start_time  = std::chrono::steady_clock::now();

   while(true)
     {
      wave_pipe::ProcessResponse response{};
      ProcessResult local{};
      wave_pipe::Status status = wave_pipe::Status::ERROR;

      {
       std::lock_guard<std::mutex> lock(g_pipe_mutex);

       HANDLE pipe = ensure_pipe_connected();
       if(pipe == INVALID_HANDLE_VALUE)
          return alglib_gpu::JobStatus::ErrorNoBackend;

       wave_pipe::FetchRequest fetch{};
       fetch.magic    = wave_pipe::MESSAGE_MAGIC;
       fetch.version  = wave_pipe::PROTOCOL_VERSION;
       fetch.cmd      = static_cast<std::uint32_t>(wave_pipe::Command::PROCESS_FETCH);
       fetch.user_tag = tag;

       if(!pipe_write(pipe, &fetch, sizeof(fetch)))
         {
          close_pipe();
          return alglib_gpu::JobStatus::ErrorBackendFailure;
         }
       FlushFileBuffers(pipe);

       if(!pipe_read(pipe, &response, sizeof(response)))
         {
          close_pipe();
          return alglib_gpu::JobStatus::ErrorBackendFailure;
         }

       status = static_cast<wave_pipe::Status>(response.status);
       local.status = status;

       if(status == wave_pipe::Status::OK)
         {
          if(response.primary_count > 0)
            {
             local.primary.resize(response.primary_count);
             const size_t bytes = static_cast<size_t>(response.primary_count) * sizeof(double);
             if(!pipe_read(pipe, local.primary.data(), bytes))
               {
                close_pipe();
                return alglib_gpu::JobStatus::ErrorBackendFailure;
               }
            }

          if(response.secondary_count > 0)
            {
             local.secondary.resize(response.secondary_count);
             const size_t bytes = static_cast<size_t>(response.secondary_count) * sizeof(double);
             if(!pipe_read(pipe, local.secondary.data(), bytes))
               {
                close_pipe();
                return alglib_gpu::JobStatus::ErrorBackendFailure;
               }
            }

          if(response.cycle_count > 0)
            {
             local.cycles.resize(response.cycle_count);
             const size_t bytes = static_cast<size_t>(response.cycle_count) * sizeof(wave_pipe::FftCyclePayload);
             if(!pipe_read(pipe, local.cycles.data(), bytes))
               {
                close_pipe();
                return alglib_gpu::JobStatus::ErrorBackendFailure;
               }
            }

          if(response.instant_count > 0)
            {
             local.instant.resize(response.instant_count);
             const size_t bytes = static_cast<size_t>(response.instant_count) * sizeof(wave_pipe::InstantMetricsPayload);
             if(!pipe_read(pipe, local.instant.data(), bytes))
               {
                close_pipe();
                return alglib_gpu::JobStatus::ErrorBackendFailure;
               }
            }

          if(response.transition_count > 0)
            {
             local.transitions.resize(response.transition_count);
             const size_t bytes = static_cast<size_t>(response.transition_count) * sizeof(wave_pipe::TransitionPayload);
             if(!pipe_read(pipe, local.transitions.data(), bytes))
               {
                close_pipe();
                return alglib_gpu::JobStatus::ErrorBackendFailure;
               }
            }
         }
      } // lock scope

      if(status == wave_pipe::Status::PENDING)
        {
         if(wait_ms == 0)
            return alglib_gpu::JobStatus::Pending;

         if(has_timeout &&
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_time).count() > wait_ms)
           return alglib_gpu::JobStatus::ErrorTimeout;

         Sleep(1);
         continue;
        }

      clear_pending_handle(tag);
      result = std::move(local);
      return status_from_pipe(status);
     }
  }

#endif // _WIN32

alglib_gpu::JobStatus execute_process(wave_pipe::Operation operation,
                                      const void* params,
                                      std::uint32_t param_size,
                                      const double* primary,
                                      int primary_len,
                                      const double* secondary,
                                      int secondary_len,
                                      ProcessResult& result)
  {
#ifdef _WIN32
   std::uint64_t handle = 0;
   auto submit_status = submit_process_async(operation,
                                             primary,
                                             primary_len,
                                             secondary,
                                             secondary_len,
                                             params,
                                             param_size,
                                             handle);

   if(handle == 0)
      return submit_status;

   return fetch_process_async(handle, result, g_active_timeout_ms.load());
#else
   (void)operation;
   (void)params;
   (void)param_size;
   (void)primary;
   (void)primary_len;
   (void)secondary;
   (void)secondary_len;
   result = ProcessResult{};
   return alglib_gpu::JobStatus::ErrorNoBackend;
#endif
  }

alglib_gpu::JobStatus wait_for_fft_result(std::uint64_t handle,
                                          std::vector<double>& buffer,
                                          int expected_len,
                                          int timeout_ms)
  {
   buffer.resize(expected_len);
   const bool use_timeout = timeout_ms > 0;
   const auto deadline = use_timeout
                             ? std::chrono::steady_clock::now() + std::chrono::milliseconds(timeout_ms)
                             : std::chrono::steady_clock::time_point::max();

   while(true)
     {
      bool ready = false;
      const auto status = g_fft_dispatcher.collect(handle, buffer.data(), expected_len, ready);
      if(status == alglib_gpu::JobStatus::Pending || !ready)
        {
         if(use_timeout && std::chrono::steady_clock::now() >= deadline)
            return alglib_gpu::JobStatus::ErrorTimeout;
         std::this_thread::sleep_for(std::chrono::milliseconds(1));
         continue;
        }
      return status;
     }
  }

// Nota: caminho legado FFT_SUBMIT/FFT_FETCH_RESULT removido. Todas as FFTs usam PROCESS_SUBMIT/PROCESS_FETCH

alglib_gpu::GpuConfig to_gpu_config(const AlglibGpuConfig* cfg)
  {
   alglib_gpu::GpuConfig config{};
   config.backend         = alglib_gpu::BackendType::Auto;
   config.device_index    = 0;
   config.stream_count    = kDefaultStreamCount;
   config.min_gpu_window  = kDefaultMinWindow;
   config.job_timeout_ms  = g_active_timeout_ms.load();

   if(cfg == nullptr)
      return config;

   switch(cfg->backend)
     {
      case 1:
         config.backend = alglib_gpu::BackendType::CpuFallback;
         break;
      case 2:
         config.backend = alglib_gpu::BackendType::Cuda;
         break;
      case 3:
         config.backend = alglib_gpu::BackendType::OpenCL;
         break;
      default:
         config.backend = alglib_gpu::BackendType::Auto;
         break;
     }

   config.device_index   = cfg->device_index;
   config.stream_count   = (cfg->stream_count > 0) ? cfg->stream_count : kDefaultStreamCount;
   config.min_gpu_window = (cfg->min_gpu_window > 0) ? cfg->min_gpu_window : kDefaultMinWindow;
   config.job_timeout_ms = (cfg->job_timeout_ms > 0) ? cfg->job_timeout_ms : g_active_timeout_ms.load();

   return config;
  }

alglib_gpu::FftKind to_fft_kind(int value)
  {
   return (value == 1) ? alglib_gpu::FftKind::Real : alglib_gpu::FftKind::Complex;
  }
} // namespace

extern "C" {

ALGLIB_GPU_API
int ALGLIB_gpu_init(const AlglibGpuConfig* config)
  {
   // Mantido por compatibilidade: não inicializa caminho local.
   // Via única: serviço/pipe (PROCESS_SUBMIT/FETCH).
   (void)config;
   g_fft_ready.store(true);
   return static_cast<int>(alglib_gpu::JobStatus::Ok);
  }

ALGLIB_GPU_API
  void ALGLIB_gpu_shutdown()
  {
   std::lock_guard<std::mutex> lock(g_fft_mutex);
   g_fft_ready.store(false);
#ifdef _WIN32
   close_pipe();
#endif
  }




// (legacy ALGLIB_fast* FFT exports removidos)

// ---------------------------------------------------------------------------
// VIA A (padronizado) – wrappers ALGLIB_wave_fft_*
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_fft_real(const double* input,
                         int           input_len,
                         int           inverse,
                         double*       output,
                         int           output_capacity,
                         int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ProcessResult result;
   auto status = execute_process(inverse ? wave_pipe::Operation::FFT_REAL_INVERSE
                                         : wave_pipe::Operation::FFT_REAL_FORWARD,
                                 nullptr,
                                 0,
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   const int produced = static_cast<int>(result.primary.size());
   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   std::copy(result.primary.begin(), result.primary.end(), output);
   *output_len = produced;
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_real_submit(const double* input,
                                int           input_len,
                                int           inverse,
                                unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *handle_out = 0ULL;
   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::uint64_t tag = 0;
   auto status = submit_process_async(inverse ? wave_pipe::Operation::FFT_REAL_INVERSE
                                             : wave_pipe::Operation::FFT_REAL_FORWARD,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      nullptr,
                                      0,
                                      tag);
   *handle_out = tag;
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_real_fetch(unsigned long long handle,
                               double*            output,
                               int                output_capacity,
                               int*               output_len,
                               int                wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle), result, wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   const int produced = static_cast<int>(result.primary.size());
   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   std::copy(result.primary.begin(), result.primary.end(), output);
   *output_len = produced;
   return static_cast<int>(alglib_gpu::JobStatus::Ok);
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_complex(const double* input,
                            int           input_len,
                            int           inverse,
                            double*       output,
                            int           output_capacity,
                            int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;
   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ProcessResult result;
   auto status = execute_process(inverse ? wave_pipe::Operation::FFT_COMPLEX_INVERSE
                                         : wave_pipe::Operation::FFT_COMPLEX_FORWARD,
                                 nullptr,
                                 0,
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   const int produced = static_cast<int>(result.primary.size());
   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   std::copy(result.primary.begin(), result.primary.end(), output);
   *output_len = produced;
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_complex_submit(const double* input,
                                   int           input_len,
                                   int           inverse,
                                   unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *handle_out = 0ULL;
   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::uint64_t tag = 0;
   auto status = submit_process_async(inverse ? wave_pipe::Operation::FFT_COMPLEX_INVERSE
                                             : wave_pipe::Operation::FFT_COMPLEX_FORWARD,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      nullptr,
                                      0,
                                      tag);
   *handle_out = tag;
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_complex_fetch(unsigned long long handle,
                                  double*            output,
                                  int                output_capacity,
                                  int*               output_len,
                                  int                wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle), result, wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   const int produced = static_cast<int>(result.primary.size());
   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   std::copy(result.primary.begin(), result.primary.end(), output);
   *output_len = produced;
   return static_cast<int>(alglib_gpu::JobStatus::Ok);
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_sine(const double* input,
                         int           input_len,
                         int           inverse,
                         double*       output,
                         int           output_capacity,
                         int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;
   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::int32_t flag = inverse ? 1 : 0;
   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::FFT_SINE_TRANSFORM,
                                 &flag,
                                 sizeof(flag),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);
   const int produced = static_cast<int>(result.primary.size());
   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   std::copy(result.primary.begin(), result.primary.end(), output);
   *output_len = produced;
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_cosine(const double* input,
                           int           input_len,
                           int           inverse,
                           double*       output,
                           int           output_capacity,
                           int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;
   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::int32_t flag = inverse ? 1 : 0;
   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::FFT_COSINE_TRANSFORM,
                                 &flag,
                                 sizeof(flag),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);
   const int produced = static_cast<int>(result.primary.size());
   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   std::copy(result.primary.begin(), result.primary.end(), output);
   *output_len = produced;
   return static_cast<int>(status);
  }

// ---------------------------------------------------------------------------
// FFT - duas séries reais (A e B) – Via A
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_fft_two_real(const double* series_a,
                             int           length_a,
                             const double* series_b,
                             int           length_b,
                             double*       spectrum_a_out,
                             int           spectrum_a_capacity,
                             int*          spectrum_a_len,
                             double*       spectrum_b_out,
                             int           spectrum_b_capacity,
                             int*          spectrum_b_len)
  {
   if(spectrum_a_len) *spectrum_a_len = 0;
   if(spectrum_b_len) *spectrum_b_len = 0;
   if(series_a == nullptr || series_b == nullptr || length_a <= 0 || length_b <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::FFT_TWO_REAL,
                                 nullptr,
                                 0,
                                 series_a,
                                 length_a,
                                 series_b,
                                 length_b,
                                 result);

   const int out_a_count = static_cast<int>(result.primary.size());
   const int out_b_count = static_cast<int>(result.secondary.size());
   if(spectrum_a_len) *spectrum_a_len = out_a_count;
   if(spectrum_b_len) *spectrum_b_len = out_b_count;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(out_a_count > 0)
     {
      if(spectrum_a_out == nullptr || spectrum_a_capacity < out_a_count)
         return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
      std::copy(result.primary.begin(), result.primary.end(), spectrum_a_out);
     }
   if(out_b_count > 0)
     {
      if(spectrum_b_out == nullptr || spectrum_b_capacity < out_b_count)
         return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
      std::copy(result.secondary.begin(), result.secondary.end(), spectrum_b_out);
     }

   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_two_real_submit(const double* series_a,
                                    int           length_a,
                                    const double* series_b,
                                    int           length_b,
                                    unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *handle_out = 0ULL;
   if(series_a == nullptr || series_b == nullptr || length_a <= 0 || length_b <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::FFT_TWO_REAL,
                                      series_a,
                                      length_a,
                                      series_b,
                                      length_b,
                                      nullptr,
                                      0,
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_fft_two_real_fetch(unsigned long long handle,
                                   double*            spectrum_a_out,
                                   int                spectrum_a_capacity,
                                   int*               spectrum_a_len,
                                   double*            spectrum_b_out,
                                   int                spectrum_b_capacity,
                                   int*               spectrum_b_len,
                                   int                wait_ms)
  {
   if(spectrum_a_len) *spectrum_a_len = 0;
   if(spectrum_b_len) *spectrum_b_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle), result, wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto st_a = copy_primary_output(result, spectrum_a_out, spectrum_a_capacity, spectrum_a_len);
   if(st_a != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(st_a);
   auto st_b = copy_secondary_output(result, spectrum_b_out, spectrum_b_capacity, spectrum_b_len);
   return static_cast<int>(st_b);
#else
   (void)handle; (void)spectrum_a_out; (void)spectrum_a_capacity; (void)spectrum_a_len;
   (void)spectrum_b_out; (void)spectrum_b_capacity; (void)spectrum_b_len; (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

// ---------------------------------------------------------------------------
// Wave operations via GPU service
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_resample(const double* input,
                         int           input_len,
                         double        factor,
                         double        cutoff,
                         int           method,
                         double*       output,
                         int           output_capacity,
                         int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ResampleParams params{};
   params.factor   = factor;
   params.cutoff   = cutoff;
   params.method   = method;
   params.reserved = 0;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_RESAMPLE,
                                 &params,
                                 sizeof(params),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);

   const int produced = static_cast<int>(result.primary.size());
   *output_len = produced;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::copy(result.primary.begin(), result.primary.end(), output);
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_resample_submit(const double* input,
                                int           input_len,
                                double        factor,
                                double        cutoff,
                                int           method,
                                unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ResampleParams params{};
   params.factor   = factor;
   params.cutoff   = cutoff;
   params.method   = method;
   params.reserved = 0;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_RESAMPLE,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_resample_fetch(unsigned long long handle,
                               double*       output,
                               int           output_capacity,
                               int*          output_len,
                               int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_zero_pad(const double* input,
                         int           input_len,
                         int           pad_left,
                         int           pad_right,
                         double*       output,
                         int           output_capacity,
                         int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ZeroPadParams params{};
   params.pad_left  = pad_left;
   params.pad_right = pad_right;
   params.reserved1 = 0;
   params.reserved2 = 0;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_ZERO_PAD,
                                 &params,
                                 sizeof(params),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);

   const int produced = static_cast<int>(result.primary.size());
   *output_len = produced;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::copy(result.primary.begin(), result.primary.end(), output);
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_zero_pad_submit(const double* input,
                                int           input_len,
                                int           pad_left,
                                int           pad_right,
                                unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ZeroPadParams params{};
   params.pad_left  = pad_left;
   params.pad_right = pad_right;
   params.reserved1 = 0;
   params.reserved2 = 0;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_ZERO_PAD,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_zero_pad_fetch(unsigned long long handle,
                               double*       output,
                               int           output_capacity,
                               int*          output_len,
                               int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_remove_dc(const double* input,
                          int           input_len,
                          int           mode,
                          double        alpha,
                          double*       output,
                          int           output_capacity,
                          int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   RemoveDcParams params{};
   params.alpha     = alpha;
   params.mode      = mode;
   params.reserved1 = 0;
   params.reserved2 = 0;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_REMOVE_DC,
                                 &params,
                                 sizeof(params),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);

   const int produced = static_cast<int>(result.primary.size());
   *output_len = produced;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::copy(result.primary.begin(), result.primary.end(), output);
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_remove_dc_submit(const double* input,
                                 int           input_len,
                                 int           mode,
                                 double        alpha,
                                 unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   RemoveDcParams params{};
   params.alpha     = alpha;
   params.mode      = mode;
   params.reserved1 = 0;
   params.reserved2 = 0;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_REMOVE_DC,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_remove_dc_fetch(unsigned long long handle,
                                double*       output,
                                int           output_capacity,
                                int*          output_len,
                                int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_filter_gaussian(const double* input,
                                int           input_len,
                                double        sigma_low,
                                double        sigma_high,
                                double        gain,
                                int           mode,
                                double*       output,
                                int           output_capacity,
                                int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   GaussianParams params{};
   params.sigma_low  = sigma_low;
   params.sigma_high = sigma_high;
   params.gain       = gain;
   params.mode       = mode;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_FILTER_GAUSSIAN,
                                 &params,
                                 sizeof(params),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);

   const int produced = static_cast<int>(result.primary.size());
   *output_len = produced;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::copy(result.primary.begin(), result.primary.end(), output);
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_filter_gaussian_submit(const double* input,
                                       int           input_len,
                                       double        sigma_low,
                                       double        sigma_high,
                                       double        gain,
                                       int           mode,
                                       unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   GaussianParams params{};
   params.sigma_low  = sigma_low;
   params.sigma_high = sigma_high;
   params.gain       = gain;
   params.mode       = mode;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_FILTER_GAUSSIAN,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_filter_gaussian_fetch(unsigned long long handle,
                                      double*       output,
                                      int           output_capacity,
                                      int*          output_len,
                                      int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_filter_notch(const double* input,
                             int           input_len,
                             double        center,
                             double        width,
                             double        depth,
                             double        gain,
                             double*       output,
                             int           output_capacity,
                             int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   NotchParams params{};
   params.center = center;
   params.width  = width;
   params.depth  = depth;
   params.gain   = gain;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_FILTER_NOTCH,
                                 &params,
                                 sizeof(params),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);

   const int produced = static_cast<int>(result.primary.size());
   *output_len = produced;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(output == nullptr || output_capacity < produced)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   std::copy(result.primary.begin(), result.primary.end(), output);
   return static_cast<int>(status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_filter_notch_submit(const double* input,
                                    int           input_len,
                                    double        center,
                                    double        width,
                                    double        depth,
                                    double        gain,
                                    unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   NotchParams params{};
   params.center = center;
   params.width  = width;
   params.depth  = depth;
   params.gain   = gain;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_FILTER_NOTCH,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_filter_notch_fetch(unsigned long long handle,
                                    double*       output,
                                    int           output_capacity,
                                    int*          output_len,
                                    int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

// ---------------------------------------------------------------------------
// Apply Mask
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_apply_mask(const double* spectrum,
                           int           spectrum_len,
                           const double* mask,
                           int           mask_len,
                           int           mode,
                           double*       output,
                           int           output_capacity,
                           int*          output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(spectrum == nullptr || mask == nullptr || spectrum_len <= 0 || mask_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   MaskParams params{};
   params.mode = mode;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_APPLY_MASK,
                                 &params,
                                 sizeof(params),
                                 spectrum,
                                 spectrum_len,
                                 mask,
                                 mask_len,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_apply_mask_submit(const double* spectrum,
                                  int           spectrum_len,
                                  const double* mask,
                                  int           mask_len,
                                  int           mode,
                                  unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(spectrum == nullptr || mask == nullptr || spectrum_len <= 0 || mask_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   MaskParams params{};
   params.mode = mode;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_APPLY_MASK,
                                      spectrum,
                                      spectrum_len,
                                      mask,
                                      mask_len,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_apply_mask_fetch(unsigned long long handle,
                                 double*       output,
                                 int           output_capacity,
                                 int*          output_len,
                                 int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

// ---------------------------------------------------------------------------
// Spectral Denoise
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_denoise(const double* spectrum,
                         int          spectrum_len,
                         int          method,
                         double       threshold,
                         double       beta,
                         int          iterations,
                         double*      output,
                         int          output_capacity,
                         int*         output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(spectrum == nullptr || spectrum_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   DenoiseParams params{};
   params.method     = method;
   params.threshold  = threshold;
   params.beta       = beta;
   params.iterations = iterations;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_DENOISE,
                                 &params,
                                 sizeof(params),
                                 spectrum,
                                 spectrum_len,
                                 nullptr,
                                 0,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_denoise_submit(const double* spectrum,
                               int           spectrum_len,
                               int           method,
                               double        threshold,
                               double        beta,
                               int           iterations,
                               unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(spectrum == nullptr || spectrum_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   DenoiseParams params{};
   params.method     = method;
   params.threshold  = threshold;
   params.beta       = beta;
   params.iterations = iterations;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_DENOISE,
                                      spectrum,
                                      spectrum_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_denoise_fetch(unsigned long long handle,
                              double*       output,
                              int           output_capacity,
                              int*          output_len,
                              int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

// ---------------------------------------------------------------------------
// Spectral Upscale / Downscale
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_upscale(const double* spectrum,
                        int          spectrum_len,
                        double       factor,
                        int          mode,
                        int          normalize,
                        double*      output,
                        int          output_capacity,
                        int*         output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(spectrum == nullptr || spectrum_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   UpscaleParams params{};
   params.factor    = factor;
   params.mode      = mode;
   params.normalize = normalize;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_UPSCALE,
                                 &params,
                                 sizeof(params),
                                 spectrum,
                                 spectrum_len,
                                 nullptr,
                                 0,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_upscale_submit(const double* spectrum,
                               int           spectrum_len,
                               double        factor,
                               int           mode,
                               int           normalize,
                               unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(spectrum == nullptr || spectrum_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   UpscaleParams params{};
   params.factor    = factor;
   params.mode      = mode;
   params.normalize = normalize;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_UPSCALE,
                                      spectrum,
                                      spectrum_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_upscale_fetch(unsigned long long handle,
                              double*       output,
                              int           output_capacity,
                              int*          output_len,
                              int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_downscale(const double* spectrum,
                          int          spectrum_len,
                          double       factor,
                          int          mode,
                          int          anti_alias,
                          double*      output,
                          int          output_capacity,
                          int*         output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(spectrum == nullptr || spectrum_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   DownscaleParams params{};
   params.factor     = factor;
   params.mode       = mode;
   params.anti_alias = anti_alias;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_DOWNSCALE,
                                 &params,
                                 sizeof(params),
                                 spectrum,
                                 spectrum_len,
                                 nullptr,
                                 0,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_downscale_submit(const double* spectrum,
                                 int           spectrum_len,
                                 double        factor,
                                 int           mode,
                                 int           anti_alias,
                                 unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(spectrum == nullptr || spectrum_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   DownscaleParams params{};
   params.factor     = factor;
   params.mode       = mode;
   params.anti_alias = anti_alias;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_DOWNSCALE,
                                      spectrum,
                                      spectrum_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_downscale_fetch(unsigned long long handle,
                                double*       output,
                                int           output_capacity,
                                int*          output_len,
                                int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

// ---------------------------------------------------------------------------
// Spectral Convolution
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_convolution(const double* spectrum,
                            int          spectrum_len,
                            const double* kernel,
                            int          kernel_len,
                            int          normalize,
                            double*      output,
                            int          output_capacity,
                            int*         output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(spectrum == nullptr || kernel == nullptr || spectrum_len <= 0 || kernel_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ConvolutionParams params{};
   params.normalize = normalize;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_CONVOLUTION,
                                 &params,
                                 sizeof(params),
                                 spectrum,
                                 spectrum_len,
                                 kernel,
                                 kernel_len,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_convolution_submit(const double* spectrum,
                                   int           spectrum_len,
                                   const double* kernel,
                                   int           kernel_len,
                                   int           normalize,
                                   unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(spectrum == nullptr || kernel == nullptr || spectrum_len <= 0 || kernel_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ConvolutionParams params{};
   params.normalize = normalize;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_CONVOLUTION,
                                      spectrum,
                                      spectrum_len,
                                      kernel,
                                      kernel_len,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_convolution_fetch(unsigned long long handle,
                                  double*       output,
                                  int           output_capacity,
                                  int*          output_len,
                                  int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

// ---------------------------------------------------------------------------
// Spectral Correlation
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_correlation(const double* spectrum,
                            int          spectrum_len,
                            const double* pattern,
                            int          pattern_len,
                            double*      output,
                            int          output_capacity,
                            int*         output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(spectrum == nullptr || pattern == nullptr || spectrum_len <= 0 || pattern_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_CORRELATION,
                                 nullptr,
                                 0,
                                 spectrum,
                                 spectrum_len,
                                 pattern,
                                 pattern_len,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_correlation_submit(const double* spectrum,
                                   int           spectrum_len,
                                   const double* pattern,
                                   int           pattern_len,
                                   unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(spectrum == nullptr || pattern == nullptr || spectrum_len <= 0 || pattern_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_CORRELATION,
                                      spectrum,
                                      spectrum_len,
                                      pattern,
                                      pattern_len,
                                      nullptr,
                                      0,
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)spectrum;
   (void)spectrum_len;
   (void)pattern;
   (void)pattern_len;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_correlation_fetch(unsigned long long handle,
                                  double*       output,
                                  int           output_capacity,
                                  int*          output_len,
                                  int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

// ---------------------------------------------------------------------------
// Phase Unwrap
// ---------------------------------------------------------------------------
ALGLIB_GPU_API
int ALGLIB_wave_phase_unwrap(const double* spectrum,
                             int          spectrum_len,
                             int          method,
                             double*      output,
                             int          output_capacity,
                             int*         output_len)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

   if(spectrum == nullptr || spectrum_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   PhaseParams params{};
   params.method = method;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_PHASE_UNWRAP,
                                 &params,
                                 sizeof(params),
                                 spectrum,
                                 spectrum_len,
                                 nullptr,
                                 0,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
  }

ALGLIB_GPU_API
int ALGLIB_wave_phase_unwrap_submit(const double* spectrum,
                                    int           spectrum_len,
                                    int           method,
                                    unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(spectrum == nullptr || spectrum_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   PhaseParams params{};
   params.method = method;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_PHASE_UNWRAP,
                                      spectrum,
                                      spectrum_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_phase_unwrap_fetch(unsigned long long handle,
                                   double*       output,
                                   int           output_capacity,
                                   int*          output_len,
                                   int           wait_ms)
  {
   if(output_len == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
   *output_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_primary_output(result, output, output_capacity, output_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)output;
   (void)output_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_analyze_spectrum(const double* input,
                                 int           input_len,
                                 int           top_k,
                                 double        min_period,
                                 double        max_period,
                                 double        energy_threshold,
                                 double        upscale_factor,
                                 double*       spectrum_out,
                                 int           spectrum_capacity,
                                 int*          spectrum_len,
                                 double*       cycles_out,
                                 int           cycles_stride,
                                 int           cycles_capacity,
                                 int*          cycles_len,
                                 double*       instant_out,
                                 int           instant_stride,
                                 int           instant_capacity,
                                 int*          instant_len)
  {
   if(spectrum_len) *spectrum_len = 0;
   if(cycles_len)   *cycles_len   = 0;
   if(instant_len)  *instant_len  = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   SpectrumParams params{};
   params.top_k           = top_k;
   params.min_period      = min_period;
   params.max_period      = max_period;
   params.energy_threshold= energy_threshold;
   params.upscale_factor  = upscale_factor;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_ANALYZE,
                                 &params,
                                 sizeof(params),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);

   const int spectrum_count = static_cast<int>(result.primary.size());
   const int cycles_count   = static_cast<int>(result.cycles.size());
   const int instant_count  = static_cast<int>(result.instant.size());

   if(spectrum_len) *spectrum_len = spectrum_count;
   if(cycles_len)   *cycles_len   = cycles_count;
   if(instant_len)  *instant_len  = instant_count;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(spectrum_count > 0)
     {
      if(spectrum_out == nullptr || spectrum_capacity < spectrum_count)
         return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
      std::copy(result.primary.begin(), result.primary.end(), spectrum_out);
     }
   else if(spectrum_out != nullptr && spectrum_capacity > 0)
     {
      // nothing to copy
     }

   if(cycles_count > 0)
     {
      if(cycles_out == nullptr || cycles_stride < 5 || cycles_capacity < cycles_count)
         return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
      for(int i = 0; i < cycles_count; ++i)
        {
         const auto& cycle = result.cycles[static_cast<std::size_t>(i)];
         double* row = cycles_out + static_cast<std::size_t>(i) * cycles_stride;
         row[0] = static_cast<double>(cycle.index);
         row[1] = cycle.amplitude;
         row[2] = cycle.phase;
         row[3] = cycle.frequency;
         row[4] = cycle.period;
        }
     }

   if(instant_count > 0)
     {
      if(instant_out == nullptr || instant_stride < 6 || instant_capacity < instant_count)
         return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
      for(int i = 0; i < instant_count; ++i)
        {
         const auto& inst = result.instant[static_cast<std::size_t>(i)];
         double* row = instant_out + static_cast<std::size_t>(i) * instant_stride;
         row[0] = inst.amplitude;
         row[1] = inst.frequency;
         row[2] = inst.energy;
         row[3] = inst.phase;
         row[4] = inst.envelope_up;
         row[5] = inst.envelope_down;
        }
     }

  return static_cast<int>(status);
 }

ALGLIB_GPU_API
int ALGLIB_wave_analyze_spectrum_submit(const double* input,
                                        int           input_len,
                                        int           top_k,
                                        double        min_period,
                                        double        max_period,
                                        double        energy_threshold,
                                        double        upscale_factor,
                                        unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   SpectrumParams params{};
   params.top_k           = top_k;
   params.min_period      = min_period;
   params.max_period      = max_period;
   params.energy_threshold= energy_threshold;
   params.upscale_factor  = upscale_factor;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_ANALYZE,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_analyze_spectrum_fetch(unsigned long long handle,
                                       double*       spectrum_out,
                                       int           spectrum_capacity,
                                       int*          spectrum_len,
                                       double*       complex_out,
                                       int           complex_capacity,
                                       int*          complex_len,
                                       double*       cycles_out,
                                       int           cycles_stride,
                                       int           cycles_capacity,
                                       int*          cycles_len,
                                       double*       instant_out,
                                       int           instant_stride,
                                       int           instant_capacity,
                                       int*          instant_len,
                                       int           wait_ms)
  {
   if(spectrum_len) *spectrum_len = 0;
   if(complex_len)  *complex_len  = 0;
   if(cycles_len)   *cycles_len   = 0;
   if(instant_len)  *instant_len  = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto primary_status = copy_primary_output(result, spectrum_out, spectrum_capacity, spectrum_len);
   if(primary_status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(primary_status);

   auto complex_status = copy_secondary_output(result, complex_out, complex_capacity, complex_len);
   if(complex_status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(complex_status);

   auto cycles_status = copy_cycles_output(result, cycles_out, cycles_stride, cycles_capacity, cycles_len);
   if(cycles_status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(cycles_status);

   auto instant_status = copy_instant_output(result, instant_out, instant_stride, instant_capacity, instant_len);
   return static_cast<int>(instant_status);
#else
   (void)handle;
   (void)spectrum_out;
   (void)spectrum_capacity;
   (void)complex_out;
   (void)complex_capacity;
   (void)complex_len;
   (void)cycles_out;
   (void)cycles_stride;
   (void)cycles_capacity;
   (void)instant_out;
   (void)instant_stride;
   (void)instant_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_detect_transitions(const double* input,
                                   int           input_len,
                                   double        energy_threshold,
                                   double        phase_jump_threshold,
                                   int           min_duration,
                                   double*       transitions_out,
                                   int           transition_stride,
                                   int           transition_capacity,
                                   int*          transition_len)
  {
   if(transition_len) *transition_len = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   TransitionParams params{};
   params.energy_threshold     = energy_threshold;
   params.phase_jump_threshold = phase_jump_threshold;
   params.min_duration         = min_duration;
   params.type_mask            = -1;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_DETECT_TRANSITIONS,
                                 &params,
                                 sizeof(params),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);

   const int transition_count = static_cast<int>(result.transitions.size());
   if(transition_len) *transition_len = transition_count;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(transition_count > 0)
     {
      if(transitions_out == nullptr || transition_stride < 5 || transition_capacity < transition_count)
         return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
      for(int i = 0; i < transition_count; ++i)
        {
         const auto& tr = result.transitions[static_cast<std::size_t>(i)];
         double* row = transitions_out + static_cast<std::size_t>(i) * transition_stride;
         row[0] = static_cast<double>(tr.start_index);
         row[1] = static_cast<double>(tr.end_index);
         row[2] = tr.energy_ratio;
         row[3] = tr.phase_shift;
         row[4] = static_cast<double>(tr.type);
        }
     }

  return static_cast<int>(status);
 }

ALGLIB_GPU_API
int ALGLIB_wave_detect_transitions_submit(const double* input,
                                          int           input_len,
                                          double        energy_threshold,
                                          double        phase_jump_threshold,
                                          int           min_duration,
                                          unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   TransitionParams params{};
   params.energy_threshold     = energy_threshold;
   params.phase_jump_threshold = phase_jump_threshold;
   params.min_duration         = min_duration;
   params.type_mask            = -1;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_DETECT_TRANSITIONS,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_detect_transitions_fetch(unsigned long long handle,
                                         double*       transitions_out,
                                         int           transition_stride,
                                         int           transition_capacity,
                                         int*          transition_len,
                                         int           wait_ms)
  {
   if(transition_len) *transition_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_transitions_output(result,
                                              transitions_out,
                                              transition_stride,
                                              transition_capacity,
                                              transition_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)transitions_out;
   (void)transition_stride;
   (void)transition_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_instant_metrics(const double* input,
                                int           input_len,
                                int           smooth_window,
                                double        epsilon,
                                double*       instant_out,
                                int           instant_stride,
                                int           instant_capacity,
                                int*          instant_len)
  {
   if(instant_len) *instant_len = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   InstantParams params{};
   params.smooth_window = smooth_window;
   params.epsilon       = epsilon;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SPECTRAL_INSTANT_METRICS,
                                 &params,
                                 sizeof(params),
                                 input,
                                 input_len,
                                 nullptr,
                                 0,
                                 result);

   const int instant_count = static_cast<int>(result.instant.size());
   if(instant_len) *instant_len = instant_count;

   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   if(instant_count > 0)
     {
      if(instant_out == nullptr || instant_stride < 6 || instant_capacity < instant_count)
         return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);
      for(int i = 0; i < instant_count; ++i)
        {
         const auto& inst = result.instant[static_cast<std::size_t>(i)];
         double* row = instant_out + static_cast<std::size_t>(i) * instant_stride;
         row[0] = inst.amplitude;
         row[1] = inst.frequency;
         row[2] = inst.energy;
         row[3] = inst.phase;
         row[4] = inst.envelope_up;
         row[5] = inst.envelope_down;
        }
     }

  return static_cast<int>(status);
 }

ALGLIB_GPU_API
int ALGLIB_wave_instant_metrics_submit(const double* input,
                                       int           input_len,
                                       int           smooth_window,
                                       double        epsilon,
                                       unsigned long long* handle_out)
  {
   if(handle_out == nullptr)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   *handle_out = 0;

   if(input == nullptr || input_len <= 0)
      return static_cast<int>(alglib_gpu::JobStatus::ErrorInvalidArgs);

   InstantParams params{};
   params.smooth_window = smooth_window;
   params.epsilon       = epsilon;

#ifdef _WIN32
   std::uint64_t handle = 0;
   auto status = submit_process_async(wave_pipe::Operation::SPECTRAL_INSTANT_METRICS,
                                      input,
                                      input_len,
                                      nullptr,
                                      0,
                                      &params,
                                      sizeof(params),
                                      handle);
   *handle_out = handle;
   return static_cast<int>(status);
#else
   (void)params;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

ALGLIB_GPU_API
int ALGLIB_wave_instant_metrics_fetch(unsigned long long handle,
                                      double*       instant_out,
                                      int           instant_stride,
                                      int           instant_capacity,
                                      int*          instant_len,
                                      int           wait_ms)
  {
   if(instant_len) *instant_len = 0;

#ifdef _WIN32
   ProcessResult result;
   auto status = fetch_process_async(static_cast<std::uint64_t>(handle),
                                     result,
                                     wait_ms);
   if(status != alglib_gpu::JobStatus::Ok)
      return static_cast<int>(status);

   auto copy_status = copy_instant_output(result,
                                          instant_out,
                                          instant_stride,
                                          instant_capacity,
                                          instant_len);
   return static_cast<int>(copy_status);
#else
   (void)handle;
   (void)instant_out;
   (void)instant_stride;
   (void)instant_capacity;
   (void)wait_ms;
   return static_cast<int>(alglib_gpu::JobStatus::ErrorNoBackend);
#endif
  }

// (legacy ALGLIB_fastsinetransform removido)

// (legacy ALGLIB_fastcosinetransform removido)

// (legacy ALGLIB_tworealffts removido)

// (legacy ALGLIB_fastconvolution removido)

// (legacy ALGLIB_fastcorellation removido)

// ======================================================================
// SSA via serviço
// ======================================================================
ALGLIB_GPU_API ssamodel* ALGLIB_ssa_create()
  {
   SsaCreateParams params{};
   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SSA_CREATE,
                                 &params,
                                 sizeof(params),
                                 nullptr,
                                 0,
                                 nullptr,
                                 0,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok || result.primary.empty())
      return nullptr;

   const std::uint64_t handle = decode_handle(result.primary[0]);
   if(handle == 0)
      return nullptr;

   ssamodel* model = new ssamodel;
   register_ssa_handle(model, handle);
   return model;
  }

ALGLIB_GPU_API void ALGLIB_ssa_destroy(ssamodel* model)
  {
   if(model == nullptr)
      return;

   std::uint64_t handle = 0;
   if(resolve_ssa_handle(model, handle) && handle != 0)
     {
      SsaHandleParam params{};
      params.handle = static_cast<double>(handle);
      ProcessResult result;
      (void)execute_process(wave_pipe::Operation::SSA_DESTROY,
                            &params,
                            sizeof(params),
                            nullptr,
                            0,
                            nullptr,
                            0,
                            result);
     }

   unregister_ssa_handle(model);
   delete model;
  }

ALGLIB_GPU_API void ALGLIB_ssa_set_window(ssamodel* model, int windowWidth)
  {
   if(model == nullptr || windowWidth <= 0)
      return;

   std::uint64_t handle = 0;
   if(!resolve_ssa_handle(model, handle) || handle == 0)
      return;

   SsaConfigParams params{};
   params.handle = static_cast<double>(handle);
   params.window = windowWidth;

   ProcessResult result;
   (void)execute_process(wave_pipe::Operation::SSA_CONFIGURE,
                         &params,
                         sizeof(params),
                         nullptr,
                         0,
                         nullptr,
                         0,
                         result);
  }

ALGLIB_GPU_API void ALGLIB_ssa_set_algo_topk_direct(ssamodel* model, int topK)
  {
   if(model == nullptr || topK <= 0)
      return;

   std::uint64_t handle = 0;
   if(!resolve_ssa_handle(model, handle) || handle == 0)
      return;

   SsaConfigParams params{};
   params.handle = static_cast<double>(handle);
   params.topk_direct = topK;

   ProcessResult result;
   (void)execute_process(wave_pipe::Operation::SSA_CONFIGURE,
                         &params,
                         sizeof(params),
                         nullptr,
                         0,
                         nullptr,
                         0,
                         result);
  }

ALGLIB_GPU_API void ALGLIB_ssa_set_algo_topk_realtime(ssamodel* model, int topK)
  {
   if(model == nullptr || topK <= 0)
      return;

   std::uint64_t handle = 0;
   if(!resolve_ssa_handle(model, handle) || handle == 0)
      return;

   SsaConfigParams params{};
   params.handle = static_cast<double>(handle);
   params.topk_realtime = topK;

   ProcessResult result;
   (void)execute_process(wave_pipe::Operation::SSA_CONFIGURE,
                         &params,
                         sizeof(params),
                         nullptr,
                         0,
                         nullptr,
                         0,
                         result);
  }

ALGLIB_GPU_API void ALGLIB_ssa_clear(ssamodel* model)
  {
   if(model == nullptr)
      return;

   std::uint64_t handle = 0;
   if(!resolve_ssa_handle(model, handle) || handle == 0)
      return;

   SsaHandleParam params{};
   params.handle = static_cast<double>(handle);

   ProcessResult result;
   (void)execute_process(wave_pipe::Operation::SSA_CLEAR,
                         &params,
                         sizeof(params),
                         nullptr,
                         0,
                         nullptr,
                         0,
                         result);
  }

ALGLIB_GPU_API
void ALGLIB_ssa_add_sequence(ssamodel* model, const double* data, int n)
  {
   if(model == nullptr || data == nullptr || n <= 0)
      return;

    std::uint64_t handle = 0;
    if(!resolve_ssa_handle(model, handle) || handle == 0)
       return;

    SsaHandleParam params{};
    params.handle = static_cast<double>(handle);

    ProcessResult result;
    (void)execute_process(wave_pipe::Operation::SSA_ADD_SEQUENCE,
                          &params,
                          sizeof(params),
                          data,
                          n,
                          nullptr,
                          0,
                          result);
  }

ALGLIB_GPU_API
void ALGLIB_ssa_append_point_and_update(ssamodel* model, double val, double updateIts)
  {
   if(model == nullptr)
      return;

   std::uint64_t handle = 0;
   if(!resolve_ssa_handle(model, handle) || handle == 0)
      return;

   SsaAppendParams params{};
   params.handle = static_cast<double>(handle);
   params.value = val;
   params.update_iters = updateIts;

   ProcessResult result;
   (void)execute_process(wave_pipe::Operation::SSA_APPEND_POINT,
                         &params,
                         sizeof(params),
                         nullptr,
                         0,
                         nullptr,
                         0,
                         result);
  }

ALGLIB_GPU_API
void ALGLIB_ssa_analyze_last(ssamodel* model, int nticks,
                             double* trendOut, double* noiseOut)
  {
   if(model == nullptr || nticks <= 0 || trendOut == nullptr || noiseOut == nullptr)
      return;

   std::uint64_t handle = 0;
   if(!resolve_ssa_handle(model, handle) || handle == 0)
      return;

   SsaAnalyzeParams params{};
   params.handle = static_cast<double>(handle);
   params.nticks = nticks;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SSA_ANALYZE_LAST,
                                 &params,
                                 sizeof(params),
                                 nullptr,
                                 0,
                                 nullptr,
                                 0,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return;

   const int trend_count = std::min(nticks, static_cast<int>(result.primary.size()));
   for(int i = 0; i < trend_count; ++i)
      trendOut[i] = result.primary[static_cast<std::size_t>(i)];

   const int noise_count = std::min(nticks, static_cast<int>(result.secondary.size()));
   for(int i = 0; i < noise_count; ++i)
      noiseOut[i] = result.secondary[static_cast<std::size_t>(i)];
  }

ALGLIB_GPU_API
void ALGLIB_ssa_forecast_last(ssamodel* model, int nticks, double* forecastOut)
  {
   if(model == nullptr || nticks <= 0 || forecastOut == nullptr)
      return;

   std::uint64_t handle = 0;
   if(!resolve_ssa_handle(model, handle) || handle == 0)
      return;

   SsaForecastParams params{};
   params.handle = static_cast<double>(handle);
   params.nticks = nticks;

   ProcessResult result;
   auto status = execute_process(wave_pipe::Operation::SSA_FORECAST_LAST,
                                 &params,
                                 sizeof(params),
                                 nullptr,
                                 0,
                                 nullptr,
                                 0,
                                 result);

   if(status != alglib_gpu::JobStatus::Ok)
      return;

   const int forecast_count = std::min(nticks, static_cast<int>(result.primary.size()));
   for(int i = 0; i < forecast_count; ++i)
      forecastOut[i] = result.primary[static_cast<std::size_t>(i)];
  }

} // extern "C"
