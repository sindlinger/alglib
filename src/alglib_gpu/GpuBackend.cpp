#include "GpuBackend.h"

#include "fasttransforms.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace alglib;

namespace alglib_gpu
{
BackendContext::BackendContext()
    : resolved_backend_(BackendType::CpuFallback),
      initialized_(false)
  {
  }

BackendContext::~BackendContext()
  {
   shutdown();
  }

JobStatus BackendContext::initialize(const GpuConfig& cfg)
  {
   std::lock_guard<std::mutex> lock(mutex_);

   if(initialized_)
      return JobStatus::Ok;

   config_           = cfg;
   resolved_backend_ = cfg.backend;

#ifdef ALGLIB_GPU_ENABLE_CUDA
   if(cfg.backend == BackendType::Cuda || cfg.backend == BackendType::Auto)
     {
      executor_ = CreateCudaExecutor(cfg);
      kernel_executor_ = CreateCudaKernelExecutor(cfg);
      if(executor_)
         resolved_backend_ = BackendType::Cuda;
      else if(cfg.backend == BackendType::Cuda)
         return JobStatus::ErrorBackendFailure;
     }
#else
   (void)cfg;
#endif

   if((!executor_ || !kernel_executor_) && (cfg.backend == BackendType::OpenCL || cfg.backend == BackendType::Auto))
     {
      if(!executor_)
         executor_ = CreateOpenClExecutor(cfg);
      if(!kernel_executor_)
         kernel_executor_ = CreateOpenClKernelExecutor(cfg);
      if(executor_ && kernel_executor_)
         resolved_backend_ = BackendType::OpenCL;
      else if(cfg.backend == BackendType::OpenCL)
         return JobStatus::ErrorBackendFailure;
     }

   if(!executor_ || !kernel_executor_)
      return JobStatus::ErrorNoBackend;

   initialized_ = true;
   return JobStatus::Ok;
  }

void BackendContext::shutdown()
  {
   std::lock_guard<std::mutex> lock(mutex_);
   executor_.reset();
   kernel_executor_.reset();
   initialized_ = false;
  }

JobStatus BackendContext::execute_fft(const FftSubmitDesc& desc, std::vector<double>& buffer)
  {
   IFftExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !executor_)
       return JobStatus::ErrorNoBackend;
    executor = executor_.get();
   }

   if(desc.kind == FftKind::Complex)
     {
      buffer.assign(desc.data, desc.data + 2 * desc.fft_len);
      return executor->execute_complex(buffer, desc.fft_len, desc.inverse);
     }

   const int expected = desc.inverse ? desc.fft_len * 2 : desc.fft_len;
   std::vector<double> input(desc.data, desc.data + expected);
   return executor->execute_real(input, buffer, desc.fft_len, desc.inverse);
 }

JobStatus BackendContext::apply_mask(const double* spectrum,
                                     std::size_t   bin_count,
                                     const double* mask,
                                     std::size_t   mask_len,
                                     bool          mask_is_complex,
                                     int           mode,
                                     std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->apply_mask(spectrum, bin_count, mask, mask_len, mask_is_complex, mode, out);
  }

JobStatus BackendContext::spectral_denoise(const double* spectrum,
                                           std::size_t   bin_count,
                                           int           method,
                                           double        threshold,
                                           double        beta,
                                           int           iterations,
                                           std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->denoise(spectrum, bin_count, method, threshold, beta, iterations, out);
  }

JobStatus BackendContext::spectral_upscale(const double* spectrum,
                                           std::size_t   bin_count,
                                           double        factor,
                                           int           mode,
                                           int           normalize,
                                           std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->upscale(spectrum, bin_count, factor, mode, normalize, out);
  }

JobStatus BackendContext::spectral_downscale(const double* spectrum,
                                             std::size_t   bin_count,
                                             double        factor,
                                             int           mode,
                                             int           anti_alias,
                                             std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->downscale(spectrum, bin_count, factor, mode, anti_alias, out);
  }

JobStatus BackendContext::spectral_convolution(const double* lhs,
                                               const double* rhs,
                                               std::size_t   bin_count,
                                               int           normalize,
                                               std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->convolution(lhs, rhs, bin_count, normalize, out);
  }

JobStatus BackendContext::spectral_correlation(const double* lhs,
                                               const double* rhs,
                                               std::size_t   bin_count,
                                               std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->correlation(lhs, rhs, bin_count, out);
  }

JobStatus BackendContext::resample_time_series(const double* input,
                                               std::size_t   length,
                                               double        factor,
                                               double        cutoff,
                                               int           method,
                                               std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->resample(input, length, factor, cutoff, method, out);
  }

JobStatus BackendContext::zero_pad_time_series(const double* input,
                                               std::size_t   length,
                                               std::size_t   pad_left,
                                               std::size_t   pad_right,
                                               std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->zero_pad(input, length, pad_left, pad_right, out);
  }

JobStatus BackendContext::remove_dc_time_series(const double* input,
                                                std::size_t   length,
                                                int           mode,
                                                double        alpha,
                                                std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->remove_dc(input, length, mode, alpha, out);
  }

JobStatus BackendContext::spectral_phase_unwrap(const double* spectrum,
                                                std::size_t   bin_count,
                                                int           method,
                                                std::vector<double>& out)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->phase_unwrap(spectrum, bin_count, method, out);
  }

JobStatus BackendContext::spectral_instant_metrics(const double* time_series,
                                                   std::size_t   length,
                                                   int           smooth_window,
                                                   double        epsilon,
                                                   std::vector<wave_pipe::InstantMetricsPayload>& instant)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->instant_metrics(time_series, length, smooth_window, epsilon, instant);
  }

JobStatus BackendContext::spectral_detect_transitions(const double* metrics,
                                                      std::size_t   length,
                                                      double        energy_threshold,
                                                      double        phase_jump_threshold,
                                                      int           min_duration,
                                                      int           type_mask,
                                                      std::vector<wave_pipe::TransitionPayload>& transitions)
  {
   IKernelExecutor* executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    executor = kernel_executor_.get();
   }

   return executor->detect_transitions(metrics,
                                       length,
                                       energy_threshold,
                                       phase_jump_threshold,
                                       min_duration,
                                       type_mask,
                                       transitions);
 }

JobStatus BackendContext::hartley_transform(const double* input,
                                            std::size_t   length,
                                            bool          inverse,
                                            std::vector<double>& out)
  {
   if(length == 0)
     {
      out.clear();
      return JobStatus::Ok;
     }

   IFftExecutor*    fft_executor    = nullptr;
   IKernelExecutor* kernel_executor = nullptr;
   {
    std::lock_guard<std::mutex> lock(mutex_);
    if(!initialized_ || !executor_ || !kernel_executor_)
       return JobStatus::ErrorNoBackend;
    fft_executor    = executor_.get();
    kernel_executor = kernel_executor_.get();
   }

   std::vector<double> real_input(input, input + length);
   std::vector<double> complex_data;
   auto status = fft_executor->execute_real(real_input,
                                            complex_data,
                                            static_cast<int>(length),
                                            false);
   if(status != JobStatus::Ok)
      return status;

   const double scale = inverse ? (1.0 / static_cast<double>(length)) : 1.0;
   return kernel_executor->hartley_from_complex(complex_data.data(),
                                                length,
                                                scale,
                                                out);
  }

BackendType BackendContext::resolved_backend() const
  {
   return resolved_backend_;
  }
} // namespace alglib_gpu

#ifndef ALGLIB_GPU_ENABLE_CUDA
std::unique_ptr<alglib_gpu::IFftExecutor> alglib_gpu::CreateCudaExecutor(const GpuConfig&)
  {
   return nullptr;
  }

std::unique_ptr<alglib_gpu::IKernelExecutor> alglib_gpu::CreateCudaKernelExecutor(const GpuConfig&)
  {
   return nullptr;
  }
#endif

#ifndef ALGLIB_GPU_ENABLE_OPENCL
std::unique_ptr<alglib_gpu::IKernelExecutor> alglib_gpu::CreateOpenClKernelExecutor(const GpuConfig&)
  {
   return nullptr;
  }
#endif
