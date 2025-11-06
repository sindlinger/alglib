#include "alglib_service_runtime.h"

#include "ap.h"
#include "dataanalysis.h"
#include "fasttransforms.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <string>
#include <ctime>
#include <cstdio>

namespace
{
using gpu_service::ProcessResultPayload;
// Legacy FFT payload removido

constexpr double kPi = 3.141592653589793238462643383279502884;

std::string escape_json(const std::string& value)
  {
   std::string escaped;
   escaped.reserve(value.size() + 8);
   for(char c : value)
     {
      switch(c)
        {
         case '\\': escaped += "\\\\"; break;
         case '\"': escaped += "\\\""; break;
         case '\b': escaped += "\\b"; break;
         case '\f': escaped += "\\f"; break;
         case '\n': escaped += "\\n"; break;
         case '\r': escaped += "\\r"; break;
         case '\t': escaped += "\\t"; break;
         default:
           if(static_cast<unsigned char>(c) < 0x20)
             {
              char buf[7];
              std::snprintf(buf, sizeof(buf), "\\u%04X", static_cast<unsigned char>(c));
              escaped += buf;
             }
           else
              escaped += c;
           break;
        }
     }
   return escaped;
  }

std::string current_timestamp_iso8601()
  {
   std::time_t now = std::time(nullptr);
   std::tm tm{};
#ifdef _WIN32
   gmtime_s(&tm, &now);
#else
   gmtime_r(&now, &tm);
#endif
   char buffer[32];
   if(std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &tm) == 0)
      return {};
   return std::string(buffer);
  }

std::string build_metadata_json(const gpu_service::ServiceRuntime::Metadata& meta)
  {
   std::ostringstream oss;
   oss << '{';
   bool first = true;
   auto append_field = [&oss, &first](const char* name, const std::string& value)
     {
      if(value.empty())
         return;
      if(!first)
         oss << ',';
      first = false;
      oss << '"' << name << '"' << ':' << '"' << escape_json(value) << '"';
     };

   append_field("version", meta.version);
   append_field("backend", meta.backend);
   append_field("driver", meta.driver);
   append_field("build", meta.build_timestamp);
   append_field("host", meta.host);

   if(meta.workers > 0)
     {
      if(!first) oss << ',';
      first = false;
      oss << "\"workers\":" << meta.workers;
     }

   if(meta.pid != 0)
     {
      if(!first) oss << ',';
      first = false;
      oss << "\"pid\":" << meta.pid;
     }

   const std::string timestamp = current_timestamp_iso8601();
   if(!timestamp.empty())
     {
      if(!first) oss << ',';
      first = false;
      oss << '"' << "timestamp" << '"' << ':' << '"' << escape_json(timestamp) << '"';
     }

   oss << '}';
   return oss.str();
  }

struct ResampleParams
  {
   double factor    = 1.0;
   double cutoff    = 0.0;
   int    method    = 0;
   int    reserved  = 0;
  };

struct ZeroPadParams
  {
   int pad_left  = 0;
   int pad_right = 0;
   int reserved1 = 0;
   int reserved2 = 0;
  };

struct RemoveDcParams
  {
   double alpha = 0.0;
   int    mode  = 0;
   int    reserved1 = 0;
   int    reserved2 = 0;
  };

struct GaussianParams
  {
   double sigma_low  = 0.05;
   double sigma_high = 0.15;
   double gain       = 1.0;
   int    mode       = 0; // 0 LP, 1 HP, 2 BP
  };

struct NotchParams
  {
   double center = 0.1;
   double width  = 0.02;
   double depth  = 1.0;
   double gain   = 1.0;
  };

struct MaskParams
  {
   std::int32_t mode      = 0;
   std::int32_t reserved1 = 0;
   std::int32_t reserved2 = 0;
   std::int32_t reserved3 = 0;
  };

struct SpectrumParams
  {
   int    top_k           = 16;
   double min_period      = 8.0;
   double max_period      = 2048.0;
   double energy_threshold= 0.0;
   double upscale_factor  = 2.0;
  };

struct TransitionParams
  {
   double energy_threshold     = 0.05;
   double phase_jump_threshold = 0.75;
   int    min_duration         = 5;
   int    type_mask            = -1;
  };

struct InstantParams
  {
   int    smooth_window = 8;
   double epsilon       = 1e-9;
  };

struct DenoiseParams
  {
   std::int32_t method     = 0;
   double       threshold  = 0.0;
   double       beta       = 0.0;
   std::int32_t iterations = 1;
  };

struct UpscaleParams
  {
   double       factor    = 1.0;
   std::int32_t mode      = 0;
   std::int32_t normalize = 0;
   std::int32_t reserved  = 0;
  };

struct DownscaleParams
  {
   double       factor     = 1.0;
   std::int32_t mode       = 0;
   std::int32_t anti_alias = 0;
   std::int32_t reserved   = 0;
  };

struct ConvolutionParams
  {
   std::int32_t normalize = 0;
   std::int32_t reserved1 = 0;
   std::int32_t reserved2 = 0;
   std::int32_t reserved3 = 0;
  };

struct PhaseParams
  {
   std::int32_t method    = 0;
   std::int32_t reserved1 = 0;
   std::int32_t reserved2 = 0;
   std::int32_t reserved3 = 0;
  };

struct SsaCreateParams
  {
   std::int32_t window        = 0;
   std::int32_t topk_direct   = 0;
   std::int32_t topk_realtime = 0;
   std::int32_t reserved      = 0;
  };

struct SsaHandleParam
  {
   double handle      = 0.0;
   double value       = 0.0;
   double extra1      = 0.0;
   double extra2      = 0.0;
  };

struct SsaConfigParams
  {
   double       handle       = 0.0;
   std::int32_t window       = 0;
   std::int32_t topk_direct  = 0;
   std::int32_t topk_realtime= 0;
  };

struct SsaAppendParams
  {
   double handle       = 0.0;
   double value        = 0.0;
   double update_iters = 0.0;
   double reserved     = 0.0;
  };

struct SsaAnalyzeParams
  {
   double       handle  = 0.0;
   std::int32_t nticks  = 0;
   std::int32_t reserved1 = 0;
   std::int32_t reserved2 = 0;
  };

struct SsaForecastParams
  {
   double       handle  = 0.0;
   std::int32_t nticks  = 0;
   std::int32_t reserved1 = 0;
   std::int32_t reserved2 = 0;
  };

template<typename T>
T unpack_params(const std::vector<std::uint8_t>& blob)
  {
   T params{};
   if(!blob.empty())
     {
      const std::size_t copy = std::min<std::size_t>(blob.size(), sizeof(T));
      std::memcpy(&params, blob.data(), copy);
     }
   return params;
  }

inline double clamp_positive(double value, double fallback)
  {
   return (value > 0.0) ? value : fallback;
  }

inline double gaussian_response(double freq, double sigma)
  {
   if(sigma <= 0.0)
      return 0.0;
   const double ratio = freq / sigma;
   return std::exp(-0.5 * ratio * ratio);
  }

void normalize_time_series(std::vector<double>& data)
  {
   if(data.empty())
      return;
   double max_abs = 0.0;
   for(double v : data)
      max_abs = std::max(max_abs, std::fabs(v));
   if(max_abs <= 0.0)
      return;
   for(double& v : data)
      v /= max_abs;
  }

std::uint64_t decode_handle(double value)
  {
   if(!std::isfinite(value) || value < 0.0)
      return 0;
   return static_cast<std::uint64_t>(std::llround(value));
  }

double encode_handle(std::uint64_t handle)
  {
   return static_cast<double>(handle);
  }

std::size_t complex_length(const std::vector<double>& buffer)
  {
   return buffer.size() / 2;
  }

bool ensure_even_complex(const std::vector<double>& buffer)
  {
   return (buffer.size() % 2) == 0;
  }

void get_complex(const std::vector<double>& buffer, std::size_t index, double& real, double& imag)
  {
   const std::size_t base = index * 2;
   real = buffer[base];
   imag = buffer[base + 1];
  }

void set_complex(std::vector<double>& buffer, std::size_t index, double real, double imag)
  {
   const std::size_t base = index * 2;
   buffer[base]     = real;
   buffer[base + 1] = imag;
  }

double complex_magnitude(const std::vector<double>& buffer, std::size_t index)
  {
   double real = 0.0;
   double imag = 0.0;
   get_complex(buffer, index, real, imag);
   return std::hypot(real, imag);
  }

} // namespace

namespace gpu_service
{
ServiceRuntime::ServiceRuntime()
    : running_(false)
    , next_ssa_handle_(1)
  {
  }

ServiceRuntime::~ServiceRuntime()
  {
   stop();
  }

bool ServiceRuntime::start(const alglib_gpu::GpuConfig& cfg, int worker_count)
  {
   if(running_)
      return true;

   auto status = backend_.initialize(cfg);
   if(status != alglib_gpu::JobStatus::Ok)
      return false;

   config_  = cfg;
   running_ = true;

   const int count = std::max(1, worker_count);
   workers_.reserve(count);
   for(int i = 0; i < count; ++i)
      workers_.emplace_back([this]() { worker_loop(); });

   return true;
  }

void ServiceRuntime::stop()
  {
   if(!running_)
      return;

   {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    running_ = false;
   }

   queue_cv_.notify_all();

   for(auto& worker : workers_)
     {
      if(worker.joinable())
         worker.join();
     }
   workers_.clear();

   backend_.shutdown();

   {
    std::lock_guard<std::mutex> lock(state_mutex_);
    pending_process_.clear();
    process_results_.clear();
    pending_batch_.clear();
   batch_results_.clear();
  }
}

void ServiceRuntime::set_metadata(const Metadata& meta)
  {
   std::lock_guard<std::mutex> lock(state_mutex_);
   metadata_ = meta;
  }

// submit_fft removido – usar submit_process com opcodes FFT_*

void ServiceRuntime::submit_process(const wave_pipe::ProcessRequest& header,
                                    const double* primary,
                                    const double* secondary,
                                    const std::uint8_t* params)
  {
   auto job = std::make_shared<ProcessJob>();
   job->header = header;
   if(header.window_len > 0 && primary)
      job->primary.assign(primary, primary + header.window_len);
   if(header.aux_len > 0 && secondary)
      job->secondary.assign(secondary, secondary + header.aux_len);
   if(header.param_size > 0 && params)
      job->params.assign(params, params + header.param_size);

   {
    std::lock_guard<std::mutex> lock(state_mutex_);
    pending_process_.insert(header.user_tag);
   }

   {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    JobWrapper wrapper{};
    wrapper.type    = JobType::Process;
    wrapper.process = job;
    queue_.push_back(wrapper);
   }

   queue_cv_.notify_one();
  }

void ServiceRuntime::submit_batch(const wave_pipe::BatchHeader& header,
                                  const wave_pipe::BatchTaskDescriptor* tasks,
                                  const std::uint8_t* data_blob,
                                  const std::uint8_t* params_blob)
  {
   auto job = std::make_shared<BatchJob>();
   job->header = header;
   if(header.task_count > 0 && tasks)
      job->tasks.assign(tasks, tasks + header.task_count);

   if(header.data_size > 0 && data_blob)
     {
      const std::size_t bytes = static_cast<std::size_t>(header.data_size);
      job->data_blob.resize(bytes);
      std::memcpy(job->data_blob.data(), data_blob, bytes);
     }

   if(header.param_size > 0 && params_blob)
     {
      const std::size_t bytes = static_cast<std::size_t>(header.param_size);
      job->params_blob.assign(params_blob, params_blob + bytes);
     }

   {
    std::lock_guard<std::mutex> lock(state_mutex_);
    pending_batch_.insert(header.user_tag);
   }

   {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    JobWrapper wrapper{};
    wrapper.type  = JobType::Batch;
    wrapper.batch = job;
    queue_.push_back(wrapper);
   }

   queue_cv_.notify_one();
  }

// fetch_fft removido – usar fetch_process

wave_pipe::Status ServiceRuntime::fetch_process(std::uint64_t tag, ProcessResultPayload& out)
  {
   std::lock_guard<std::mutex> lock(state_mutex_);
   auto it = process_results_.find(tag);
   if(it != process_results_.end())
     {
      out = it->second;
      process_results_.erase(it);
      return out.status;
     }

   if(pending_process_.count(tag) > 0)
      return wave_pipe::Status::PENDING;

   out.status = wave_pipe::Status::ERROR;
   return out.status;
  }

wave_pipe::Status ServiceRuntime::fetch_batch(std::uint64_t tag, BatchResultPayload& out)
  {
   std::lock_guard<std::mutex> lock(state_mutex_);
   auto it = batch_results_.find(tag);
   if(it != batch_results_.end())
     {
      out = it->second;
      batch_results_.erase(it);
      return out.status;
     }

   if(pending_batch_.count(tag) > 0)
      return wave_pipe::Status::PENDING;

   out.status = wave_pipe::Status::ERROR;
   return out.status;
  }

alglib_gpu::BackendType ServiceRuntime::resolved_backend() const
  {
   return backend_.resolved_backend();
  }

void ServiceRuntime::worker_loop()
  {
   while(true)
     {
      JobWrapper job;
      bool has_job = false;
      {
       std::unique_lock<std::mutex> lock(queue_mutex_);
       queue_cv_.wait(lock, [this]() { return !queue_.empty() || !running_; });
       if(!running_ && queue_.empty())
          return;

       job = queue_.front();
       queue_.pop_front();
       has_job = true;
      }

      if(!has_job)
         continue;

      switch(job.type)
        {
         // JobType::Fft removido
         case JobType::Process:
            if(job.process)
               process_signal(job.process);
            break;
         case JobType::Batch:
            if(job.batch)
               process_batch(job.batch);
            break;
        }
     }
  }

static wave_pipe::Status status_from_job(alglib_gpu::JobStatus status)
  {
   switch(status)
     {
      case alglib_gpu::JobStatus::Ok:
         return wave_pipe::Status::OK;
      case alglib_gpu::JobStatus::Pending:
         return wave_pipe::Status::PENDING;
      case alglib_gpu::JobStatus::ErrorTimeout:
         return wave_pipe::Status::TIMEOUT;
      case alglib_gpu::JobStatus::ErrorInvalidArgs:
         return wave_pipe::Status::INVALID;
      default:
         return wave_pipe::Status::ERROR;
     }
  }

wave_pipe::Status ServiceRuntime::convert_status(alglib_gpu::JobStatus status)
  {
   return status_from_job(status);
  }

// process_fft removido – FFT via PROCESS_*

static std::vector<double> build_frequency_response(const std::vector<double>& spectrum,
                                                    wave_pipe::Operation op,
                                                    const GaussianParams& g_params,
                                                    const NotchParams& n_params)
  {
   const std::size_t complex_len = spectrum.size() / 2;
   std::vector<double> response(complex_len, 1.0);
   const double nyquist = 0.5;
   for(std::size_t k = 0; k < complex_len; ++k)
     {
      const double freq = (complex_len == 0) ? 0.0 : (static_cast<double>(k) / static_cast<double>(complex_len - 1)) * nyquist;
      double gain = 1.0;
      if(op == wave_pipe::Operation::SPECTRAL_FILTER_GAUSSIAN)
        {
         const double sigma_low  = clamp_positive(g_params.sigma_low, 0.05);
         const double sigma_high = clamp_positive(g_params.sigma_high, sigma_low * 2.0);
         switch(g_params.mode)
           {
            case 0: // LP
              gain = gaussian_response(freq, sigma_low);
              break;
            case 1: // HP
              gain = 1.0 - gaussian_response(freq, sigma_low);
              break;
            case 2: // BP
            default:
              {
               const double low  = gaussian_response(freq, sigma_low);
               const double high = gaussian_response(freq, sigma_high);
               gain = std::max(0.0, high - low);
              }
              break;
           }
         gain *= g_params.gain;
        }
      else if(op == wave_pipe::Operation::SPECTRAL_FILTER_NOTCH)
        {
         const double distance = std::fabs(freq - n_params.center);
         const double width    = clamp_positive(n_params.width, 0.01);
         const double attenuation = std::exp(-0.5 * (distance / width) * (distance / width)) * n_params.depth;
         gain = std::max(0.0, n_params.gain * (1.0 - attenuation));
        }
      response[k] = gain;
     }
   return response;
  }

void ServiceRuntime::process_signal(const std::shared_ptr<ProcessJob>& job)
  {
  ProcessResultPayload result{};
  result.status = wave_pipe::Status::ERROR;

  const auto operation = static_cast<wave_pipe::Operation>(job->header.operation);
  const std::size_t window_len = job->primary.size();
  const bool is_ping = (operation == wave_pipe::Operation::PING);
  const bool requires_primary = !is_ping && !(
       operation == wave_pipe::Operation::SSA_CREATE ||
       operation == wave_pipe::Operation::SSA_DESTROY ||
       operation == wave_pipe::Operation::SSA_CONFIGURE ||
       operation == wave_pipe::Operation::SSA_CLEAR ||
       operation == wave_pipe::Operation::SSA_APPEND_POINT ||
       operation == wave_pipe::Operation::SSA_ANALYZE_LAST ||
       operation == wave_pipe::Operation::SSA_FORECAST_LAST);

   if(is_ping)
     {
      Metadata metadata_copy;
      {
       std::lock_guard<std::mutex> lock(state_mutex_);
       metadata_copy = metadata_;
      }

      result.status = wave_pipe::Status::OK;
      const std::string summary = build_metadata_json(metadata_copy);
      result.metadata.assign(summary.begin(), summary.end());
     }
   else if(requires_primary && window_len == 0)
     {
      result.status = wave_pipe::Status::INVALID;
     }
  else
    {
     // Temporariamente desativa SSA no serviço (CPU), até versão GPU estar disponível.
     if(operation >= wave_pipe::Operation::SSA_CREATE && operation <= wave_pipe::Operation::SSA_FORECAST_LAST)
       {
        result.status = wave_pipe::Status::INVALID;
       }
     else
       {
        switch(operation)
        {
        case wave_pipe::Operation::PING:
          result.status = wave_pipe::Status::OK;
          break;

         case wave_pipe::Operation::FFT_COMPLEX_FORWARD:
           {
            if(!ensure_even_complex(job->primary))
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            alglib_gpu::FftSubmitDesc desc{};
            desc.data    = job->primary.data();
            desc.fft_len = static_cast<int>(complex_length(job->primary));
            desc.kind    = alglib_gpu::FftKind::Complex;
            desc.inverse = false;

            std::vector<double> spectrum;
            auto status = backend_.execute_fft(desc, spectrum);
            result.status = convert_status(status);
            if(result.status == wave_pipe::Status::OK)
               result.primary = std::move(spectrum);
           }
           break;

         case wave_pipe::Operation::FFT_COMPLEX_INVERSE:
           {
            if(!ensure_even_complex(job->primary))
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            alglib_gpu::FftSubmitDesc desc{};
            desc.data    = job->primary.data();
            desc.fft_len = static_cast<int>(complex_length(job->primary));
            desc.kind    = alglib_gpu::FftKind::Complex;
            desc.inverse = true;

            std::vector<double> time_domain;
            auto status = backend_.execute_fft(desc, time_domain);
            result.status = convert_status(status);
            if(result.status == wave_pipe::Status::OK)
               result.primary = std::move(time_domain);
           }
           break;

         case wave_pipe::Operation::FFT_REAL_FORWARD:
           {
            if(job->primary.empty())
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            alglib_gpu::FftSubmitDesc desc{};
            desc.data    = job->primary.data();
            desc.fft_len = static_cast<int>(job->primary.size());
            desc.kind    = alglib_gpu::FftKind::Real;
            desc.inverse = false;

            std::vector<double> spectrum;
            auto status = backend_.execute_fft(desc, spectrum);
            result.status = convert_status(status);
            if(result.status == wave_pipe::Status::OK)
               result.primary = std::move(spectrum);
           }
           break;

         case wave_pipe::Operation::FFT_REAL_INVERSE:
           {
            if(!ensure_even_complex(job->primary))
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            alglib_gpu::FftSubmitDesc desc{};
            desc.data    = job->primary.data();
            desc.fft_len = static_cast<int>(complex_length(job->primary));
            desc.kind    = alglib_gpu::FftKind::Real;
            desc.inverse = true;

            std::vector<double> time_domain;
            auto status = backend_.execute_fft(desc, time_domain);
            result.status = convert_status(status);
            if(result.status == wave_pipe::Status::OK)
               result.primary = std::move(time_domain);
           }
           break;

         case wave_pipe::Operation::FFT_SINE_TRANSFORM:
         case wave_pipe::Operation::FFT_COSINE_TRANSFORM:
           {
            if(job->primary.empty())
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            bool inverse = false;
            if(job->params.size() >= sizeof(std::int32_t))
              {
               std::int32_t tmp = 0;
               std::memcpy(&tmp, job->params.data(), sizeof(std::int32_t));
               inverse = (tmp != 0);
              }

            std::vector<double> transform;
            const auto gpu_status = backend_.hartley_transform(job->primary.data(),
                                                               job->primary.size(),
                                                               inverse,
                                                               transform);
            result.status = convert_status(gpu_status);
            if(result.status == wave_pipe::Status::OK)
               result.primary = std::move(transform);
           }
           break;

         case wave_pipe::Operation::FFT_TWO_REAL:
           {
            if(job->secondary.empty() || job->primary.empty())
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            alglib_gpu::FftSubmitDesc desc_a{};
            desc_a.data    = job->primary.data();
            desc_a.fft_len = static_cast<int>(job->primary.size());
            desc_a.kind    = alglib_gpu::FftKind::Real;
            desc_a.inverse = false;

            std::vector<double> spectrum_a;
            auto status_a = backend_.execute_fft(desc_a, spectrum_a);

            alglib_gpu::FftSubmitDesc desc_b{};
            desc_b.data    = job->secondary.data();
            desc_b.fft_len = static_cast<int>(job->secondary.size());
            desc_b.kind    = alglib_gpu::FftKind::Real;
            desc_b.inverse = false;

            std::vector<double> spectrum_b;
            auto status_b = backend_.execute_fft(desc_b, spectrum_b);

            if(status_a != alglib_gpu::JobStatus::Ok)
              {
               result.status = convert_status(status_a);
               break;
              }
            if(status_b != alglib_gpu::JobStatus::Ok)
              {
               result.status = convert_status(status_b);
               break;
              }

            result.primary   = std::move(spectrum_a);
            result.secondary = std::move(spectrum_b);
            result.status    = wave_pipe::Status::OK;
           }
           break;
        case wave_pipe::Operation::SPECTRAL_RESAMPLE:
          {
           const auto params = unpack_params<ResampleParams>(job->params);
           const double factor = clamp_positive(params.factor, 1.0);
           if(factor == 1.0)
             {
              result.primary = job->primary;
              normalize_time_series(result.primary);
              result.status = wave_pipe::Status::OK;
              break;
             }

           std::vector<double> gpu_output;
           const auto gpu_status = backend_.resample_time_series(job->primary.data(),
                                                                 job->primary.size(),
                                                                 factor,
                                                                 params.cutoff,
                                                                 params.method,
                                                                 gpu_output);
           result.status = convert_status(gpu_status);
           if(result.status == wave_pipe::Status::OK)
             {
              result.primary = std::move(gpu_output);
              normalize_time_series(result.primary);
        }
       }
    }
          break;

        case wave_pipe::Operation::SPECTRAL_ZERO_PAD:
          {
           const auto params = unpack_params<ZeroPadParams>(job->params);
           const std::size_t left  = static_cast<std::size_t>(std::max(0, params.pad_left));
           const std::size_t right = static_cast<std::size_t>(std::max(0, params.pad_right));
           std::vector<double> gpu_output;
           const auto gpu_status = backend_.zero_pad_time_series(job->primary.data(),
                                                                  job->primary.size(),
                                                                  left,
                                                                  right,
                                                                  gpu_output);
           result.status = convert_status(gpu_status);
           if(result.status == wave_pipe::Status::OK)
              result.primary = std::move(gpu_output);
          }
          break;

        case wave_pipe::Operation::SPECTRAL_REMOVE_DC:
          {
           const auto params = unpack_params<RemoveDcParams>(job->params);
           std::vector<double> gpu_output;
           const auto gpu_status = backend_.remove_dc_time_series(job->primary.data(),
                                                                  job->primary.size(),
                                                                  params.mode,
                                                                  params.alpha,
                                                                  gpu_output);
           result.status = convert_status(gpu_status);
           if(result.status == wave_pipe::Status::OK)
              result.primary = std::move(gpu_output);
          }
          break;

         case wave_pipe::Operation::SPECTRAL_FILTER_GAUSSIAN:
         case wave_pipe::Operation::SPECTRAL_FILTER_NOTCH:
           {
            std::vector<double> working = job->primary;
            alglib_gpu::FftSubmitDesc desc{};
            desc.data    = working.data();
            desc.fft_len = static_cast<int>(working.size());
            desc.kind    = alglib_gpu::FftKind::Real;
            desc.inverse = false;

            std::vector<double> spectrum;
            auto status = backend_.execute_fft(desc, spectrum);
            if(status != alglib_gpu::JobStatus::Ok)
              {
               result.status = convert_status(status);
               break;
              }

            GaussianParams g_params;
            NotchParams    n_params;
            if(operation == wave_pipe::Operation::SPECTRAL_FILTER_GAUSSIAN)
               g_params = unpack_params<GaussianParams>(job->params);
            else
               n_params = unpack_params<NotchParams>(job->params);

            const auto response = build_frequency_response(spectrum, operation, g_params, n_params);
           const std::size_t complex_len = spectrum.size() / 2;

           std::vector<double> filtered_spectrum;
           const auto gpu_status = backend_.apply_mask(spectrum.data(),
                                                       complex_len,
                                                       response.data(),
                                                       response.size(),
                                                       false,
                                                       0,
                                                       filtered_spectrum);
           if(gpu_status == alglib_gpu::JobStatus::Ok)
             spectrum = std::move(filtered_spectrum);
           else
             {
              result.status = convert_status(gpu_status);
              break;
             }

            alglib_gpu::FftSubmitDesc inv{};
            inv.data    = spectrum.data();
            inv.fft_len = static_cast<int>(working.size());
            inv.kind    = alglib_gpu::FftKind::Real;
            inv.inverse = true;

            std::vector<double> filtered;
            status = backend_.execute_fft(inv, filtered);
            result.status = convert_status(status);
            if(result.status == wave_pipe::Status::OK)
               result.primary = std::move(filtered);
           }
           break;

         case wave_pipe::Operation::SPECTRAL_APPLY_MASK:
           {
            if(!ensure_even_complex(job->primary))
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            const auto params = unpack_params<MaskParams>(job->params);
            const std::size_t bins = complex_length(job->primary);

            if(job->secondary.empty())
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            const bool mask_is_complex = (job->secondary.size() == job->primary.size());
            if((mask_is_complex && job->secondary.size() != bins * 2) ||
               (!mask_is_complex && job->secondary.size() != bins))
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            std::vector<double> gpu_output;
            const auto gpu_status = backend_.apply_mask(job->primary.data(),
                                                        bins,
                                                        job->secondary.data(),
                                                        job->secondary.size(),
                                                        mask_is_complex,
                                                        params.mode,
                                                        gpu_output);
            result.status = convert_status(gpu_status);
            if(result.status == wave_pipe::Status::OK)
               result.primary = std::move(gpu_output);
           }
           break;

        case wave_pipe::Operation::SPECTRAL_DENOISE:
          {
           if(!ensure_even_complex(job->primary))
             {
              result.status = wave_pipe::Status::INVALID;
              break;
             }

           const auto params = unpack_params<DenoiseParams>(job->params);
           const std::size_t bins = complex_length(job->primary);

           std::vector<double> gpu_output;
           const auto gpu_status = backend_.spectral_denoise(job->primary.data(),
                                                              bins,
                                                              params.method,
                                                              params.threshold,
                                                              params.beta,
                                                              params.iterations,
                                                              gpu_output);
           result.status = convert_status(gpu_status);
           if(result.status == wave_pipe::Status::OK)
              result.primary = std::move(gpu_output);
          }
          break;

        case wave_pipe::Operation::SPECTRAL_UPSCALE:
          {
           if(!ensure_even_complex(job->primary))
             {
              result.status = wave_pipe::Status::INVALID;
              break;
             }

           const auto params = unpack_params<UpscaleParams>(job->params);
           const double factor = std::max(1.0, params.factor);
           const std::size_t bins = complex_length(job->primary);
           std::vector<double> gpu_output;
           const auto gpu_status = backend_.spectral_upscale(job->primary.data(),
                                                             bins,
                                                             factor,
                                                             params.mode,
                                                             params.normalize,
                                                             gpu_output);
           result.status = convert_status(gpu_status);
           if(result.status == wave_pipe::Status::OK)
              result.primary = std::move(gpu_output);
          }
          break;

        case wave_pipe::Operation::SPECTRAL_DOWNSCALE:
          {
           if(!ensure_even_complex(job->primary))
             {
              result.status = wave_pipe::Status::INVALID;
              break;
             }

           const auto params = unpack_params<DownscaleParams>(job->params);
           const double factor = std::max(1.0, params.factor);
           const std::size_t bins = complex_length(job->primary);
           std::vector<double> gpu_output;
           const auto gpu_status = backend_.spectral_downscale(job->primary.data(),
                                                               bins,
                                                               factor,
                                                               params.mode,
                                                               params.anti_alias,
                                                               gpu_output);
           result.status = convert_status(gpu_status);
           if(result.status == wave_pipe::Status::OK)
              result.primary = std::move(gpu_output);
          }
          break;

        case wave_pipe::Operation::SPECTRAL_CONVOLUTION:
          {
           if(!ensure_even_complex(job->primary) || !ensure_even_complex(job->secondary))
             {
              result.status = wave_pipe::Status::INVALID;
              break;
             }

           const auto params = unpack_params<ConvolutionParams>(job->params);
           if(job->primary.size() != job->secondary.size())
             {
              result.status = wave_pipe::Status::INVALID;
              break;
             }

           const std::size_t bins = complex_length(job->primary);
           std::vector<double> gpu_output;
           const auto gpu_status = backend_.spectral_convolution(job->primary.data(),
                                                                 job->secondary.data(),
                                                                 bins,
                                                                 params.normalize,
                                                                 gpu_output);
           result.status = convert_status(gpu_status);
           if(result.status == wave_pipe::Status::OK)
              result.primary = std::move(gpu_output);
          }
          break;

        case wave_pipe::Operation::SPECTRAL_CORRELATION:
          {
           if(!ensure_even_complex(job->primary) || !ensure_even_complex(job->secondary) || job->primary.size() != job->secondary.size())
             {
              result.status = wave_pipe::Status::INVALID;
              break;
             }

           const std::size_t bins = complex_length(job->primary);
           std::vector<double> gpu_output;
           const auto gpu_status = backend_.spectral_correlation(job->primary.data(),
                                                                 job->secondary.data(),
                                                                 bins,
                                                                 gpu_output);
           result.status = convert_status(gpu_status);
           if(result.status == wave_pipe::Status::OK)
              result.primary = std::move(gpu_output);
          }
          break;

         case wave_pipe::Operation::SPECTRAL_ANALYZE:
           {
            const auto params = unpack_params<SpectrumParams>(job->params);
            std::vector<double> working = job->primary;
            alglib_gpu::FftSubmitDesc desc{};
            desc.data    = working.data();
            desc.fft_len = static_cast<int>(working.size());
            desc.kind    = alglib_gpu::FftKind::Real;
            desc.inverse = false;

            std::vector<double> spectrum;
            auto status = backend_.execute_fft(desc, spectrum);
            if(status != alglib_gpu::JobStatus::Ok)
              {
               result.status = convert_status(status);
               break;
              }

            const std::size_t complex_len = spectrum.size() / 2;
            std::vector<std::pair<double, wave_pipe::FftCyclePayload>> ranked;
            ranked.reserve(complex_len);

            const double length = static_cast<double>(working.size());
            for(std::size_t k = 1; k < complex_len; ++k)
              {
               const double real = spectrum[2 * k];
               const double imag = spectrum[2 * k + 1];
               const double amp  = std::hypot(real, imag) / length;
               const double freq = static_cast<double>(k) / length;
               if(freq <= 0.0)
                  continue;
               const double period = (freq > 0.0) ? (1.0 / freq) : 0.0;
               if(period < params.min_period || period > params.max_period)
                  continue;
               if(amp < params.energy_threshold)
                  continue;
               wave_pipe::FftCyclePayload payload{};
               payload.index     = static_cast<std::int32_t>(k);
               payload.amplitude = amp;
               payload.phase     = std::atan2(imag, real);
               payload.frequency = freq;
               payload.period    = period;
               ranked.emplace_back(amp, payload);
              }

            std::sort(ranked.begin(), ranked.end(), [](const auto& lhs, const auto& rhs) {
             return lhs.first > rhs.first;
            });

            const std::size_t limit = std::min<std::size_t>(ranked.size(), static_cast<std::size_t>(std::max(1, params.top_k)));
            result.cycles.reserve(limit);
            for(std::size_t i = 0; i < limit; ++i)
               result.cycles.push_back(ranked[i].second);

            result.primary.resize(complex_len);
            for(std::size_t k = 0; k < complex_len; ++k)
              {
               const double real = spectrum[2 * k];
               const double imag = spectrum[2 * k + 1];
               result.primary[k] = std::hypot(real, imag) / length;
              }

            result.instant.resize(working.size());
            double prev_phase = 0.0;
            for(std::size_t i = 0; i < working.size(); ++i)
              {
               double amplitude = std::abs(working[i]);
               double phase = std::atan2(working[i], (i > 0 ? working[i - 1] : working[i]));
               double phase_diff = phase - prev_phase;
               prev_phase = phase;
               wave_pipe::InstantMetricsPayload metrics{};
               metrics.amplitude = amplitude;
               metrics.energy    = amplitude * amplitude;
               metrics.phase     = phase;
               metrics.frequency = phase_diff;
               metrics.envelope_up   = amplitude;
               metrics.envelope_down = -amplitude;
               result.instant[i] = metrics;
              }

            result.secondary = std::move(spectrum);
            result.status = wave_pipe::Status::OK;
           }
           break;

         case wave_pipe::Operation::SPECTRAL_DETECT_TRANSITIONS:
           {
            const auto params = unpack_params<TransitionParams>(job->params);
            double energy_threshold = params.energy_threshold;
            if(energy_threshold <= 0.0 && !job->primary.empty())
              {
               double total = 0.0;
               for(double v : job->primary)
                  total += v * v;
               energy_threshold = total / static_cast<double>(job->primary.size());
              }

            std::vector<wave_pipe::TransitionPayload> transitions;
            const auto status = backend_.spectral_detect_transitions(job->primary.data(),
                                                                     job->primary.size(),
                                                                     energy_threshold,
                                                                     params.phase_jump_threshold,
                                                                     params.min_duration,
                                                                     params.type_mask,
                                                                     transitions);
            result.status = convert_status(status);
            if(result.status == wave_pipe::Status::OK)
              {
               result.primary = job->primary;
               result.transitions = std::move(transitions);
              }
           }
           break;

         case wave_pipe::Operation::SPECTRAL_INSTANT_METRICS:
           {
            const auto params = unpack_params<InstantParams>(job->params);
            std::vector<wave_pipe::InstantMetricsPayload> instant;
            const auto status = backend_.spectral_instant_metrics(job->primary.data(),
                                                                  job->primary.size(),
                                                                  params.smooth_window,
                                                                  params.epsilon,
                                                                  instant);
            result.status = convert_status(status);
            if(result.status == wave_pipe::Status::OK)
              {
               result.primary = job->primary;
               result.instant = std::move(instant);
              }
           }
           break;

         case wave_pipe::Operation::SPECTRAL_PHASE_UNWRAP:
           {
            if(!ensure_even_complex(job->primary))
              {
               result.status = wave_pipe::Status::INVALID;
               break;
              }

            const auto params = unpack_params<PhaseParams>(job->params);
            std::vector<double> unwrapped;
            const auto status = backend_.spectral_phase_unwrap(job->primary.data(),
                                                               complex_length(job->primary),
                                                               params.method,
                                                               unwrapped);
            result.status = convert_status(status);
            if(result.status == wave_pipe::Status::OK)
               result.primary = std::move(unwrapped);
           }
           break;

         case wave_pipe::Operation::SSA_CREATE:
           {
            // SSA via CPU desativado até versão GPU. Mantém opcode por compatibilidade.
            result.status = wave_pipe::Status::INVALID;
           }
           break;

         case wave_pipe::Operation::SSA_DESTROY:
           {
            result.status = wave_pipe::Status::INVALID;
           }
           break;

         case wave_pipe::Operation::SSA_CONFIGURE:
           {
            result.status = wave_pipe::Status::INVALID;
           }
           break;

         case wave_pipe::Operation::SSA_CLEAR:
           {
            result.status = wave_pipe::Status::INVALID;
           }
           break;

         case wave_pipe::Operation::SSA_ADD_SEQUENCE:
           {
            result.status = wave_pipe::Status::INVALID;
           }
           break;

         case wave_pipe::Operation::SSA_APPEND_POINT:
           {
            result.status = wave_pipe::Status::INVALID;
           }
           break;

         case wave_pipe::Operation::SSA_ANALYZE_LAST:
           {
            result.status = wave_pipe::Status::INVALID;
           }
           break;

         case wave_pipe::Operation::SSA_FORECAST_LAST:
           {
            result.status = wave_pipe::Status::INVALID;
           }
           break;

         default:
           result.status = wave_pipe::Status::INVALID;
           break;
        }
     }

   {
    std::lock_guard<std::mutex> lock(state_mutex_);
    pending_process_.erase(job->header.user_tag);
    if(result.secondary.empty() && !job->secondary.empty())
       result.secondary = job->secondary;
    process_results_[job->header.user_tag] = std::move(result);
   }
  }

static void append_bytes(std::vector<std::uint8_t>& buffer, const void* data, std::size_t bytes)
  {
   if(bytes == 0)
      return;
   const std::size_t offset = buffer.size();
   buffer.resize(offset + bytes);
   std::memcpy(buffer.data() + offset, data, bytes);
  }

static std::vector<std::uint8_t> serialize_process_payload(const wave_pipe::ProcessRequest& header,
                                                           const ProcessResultPayload& payload)
  {
   wave_pipe::ProcessResponse response{};
   response.magic            = wave_pipe::MESSAGE_MAGIC;
   response.version          = wave_pipe::PROTOCOL_VERSION;
   response.cmd              = static_cast<std::uint32_t>(wave_pipe::Command::PROCESS_FETCH);
   response.user_tag         = header.user_tag;
   response.status           = static_cast<std::int32_t>(payload.status);
   response.primary_count    = static_cast<std::uint32_t>(payload.primary.size());
   response.secondary_count  = static_cast<std::uint32_t>(payload.secondary.size());
   response.cycle_count      = static_cast<std::uint32_t>(payload.cycles.size());
   response.instant_count    = static_cast<std::uint32_t>(payload.instant.size());
   response.transition_count = static_cast<std::uint32_t>(payload.transitions.size());
  response.payload_size     = static_cast<std::uint32_t>(payload.primary.size() * sizeof(double)
                                  + payload.secondary.size() * sizeof(double)
                                  + payload.cycles.size() * sizeof(wave_pipe::FftCyclePayload)
                                  + payload.instant.size() * sizeof(wave_pipe::InstantMetricsPayload)
                                  + payload.transitions.size() * sizeof(wave_pipe::TransitionPayload)
                                  + payload.metadata.size());

   std::vector<std::uint8_t> buffer;
   buffer.reserve(sizeof(response) + response.payload_size);
   append_bytes(buffer, &response, sizeof(response));
   append_bytes(buffer, payload.primary.data(), payload.primary.size() * sizeof(double));
   append_bytes(buffer, payload.secondary.data(), payload.secondary.size() * sizeof(double));
   append_bytes(buffer, payload.cycles.data(), payload.cycles.size() * sizeof(wave_pipe::FftCyclePayload));
   append_bytes(buffer, payload.instant.data(), payload.instant.size() * sizeof(wave_pipe::InstantMetricsPayload));
   append_bytes(buffer, payload.transitions.data(), payload.transitions.size() * sizeof(wave_pipe::TransitionPayload));
   append_bytes(buffer, payload.metadata.data(), payload.metadata.size());
  return buffer;
 }

void ServiceRuntime::process_batch(const std::shared_ptr<BatchJob>& job)
  {
   BatchResultPayload payload{};
   payload.status = wave_pipe::Status::OK;

   std::vector<std::uint8_t> concatenated;
   concatenated.reserve(1024);

   bool any_error = false;
   for(std::size_t idx = 0; idx < job->tasks.size(); ++idx)
     {
      const auto& task = job->tasks[idx];
      wave_pipe::ProcessRequest request{};
      request.magic     = wave_pipe::MESSAGE_MAGIC;
      request.version   = wave_pipe::PROTOCOL_VERSION;
      request.cmd       = static_cast<std::uint32_t>(wave_pipe::Command::PROCESS_SUBMIT);
      request.user_tag  = job->header.user_tag + idx + 1;
      request.operation = task.operation;
      request.flags     = task.flags;
      request.window_len= task.window_len;
      request.aux_len   = task.aux_len;
      request.param_size= static_cast<std::uint32_t>(task.param_offset < job->params_blob.size()
                          ? std::min<std::size_t>(job->params_blob.size() - static_cast<std::size_t>(task.param_offset), static_cast<std::size_t>(UINT32_MAX))
                          : 0);

      const double* primary = nullptr;
      const double* secondary = nullptr;
      const std::uint8_t* params = nullptr;

      if(task.data_offset < job->data_blob.size())
        {
         const std::size_t byte_offset = static_cast<std::size_t>(task.data_offset);
         primary = reinterpret_cast<const double*>(job->data_blob.data() + byte_offset);
        }
      if(task.aux_len > 0)
        {
         const std::size_t aux_offset = static_cast<std::size_t>(task.data_offset + static_cast<std::uint64_t>(task.window_len) * sizeof(double));
         if(aux_offset < job->data_blob.size())
            secondary = reinterpret_cast<const double*>(job->data_blob.data() + aux_offset);
        }
      if(task.param_offset < job->params_blob.size())
         params = job->params_blob.data() + static_cast<std::size_t>(task.param_offset);

      auto proc_job = std::make_shared<ProcessJob>();
      proc_job->header = request;
      if(primary && task.window_len > 0)
         proc_job->primary.assign(primary, primary + task.window_len);
      if(secondary && task.aux_len > 0)
         proc_job->secondary.assign(secondary, secondary + task.aux_len);
      if(params && request.param_size > 0)
         proc_job->params.assign(params, params + request.param_size);

      process_signal(proc_job);

      ProcessResultPayload local_result{};
      {
       std::lock_guard<std::mutex> lock(state_mutex_);
       auto it = process_results_.find(request.user_tag);
       if(it != process_results_.end())
         {
          local_result = it->second;
          process_results_.erase(it);
         }
       else
         {
          any_error = true;
          continue;
         }
      }

      const auto serialized = serialize_process_payload(request, local_result);
      append_bytes(concatenated, serialized.data(), serialized.size());
     }

   if(any_error)
      payload.status = wave_pipe::Status::ERROR;

   payload.buffer = std::move(concatenated);

   std::lock_guard<std::mutex> lock(state_mutex_);
   pending_batch_.erase(job->header.user_tag);
   batch_results_[job->header.user_tag] = std::move(payload);
  }

} // namespace gpu_service
