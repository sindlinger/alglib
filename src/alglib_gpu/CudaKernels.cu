#include "GpuBackend.h"
#include "alglib_pipe_messages.h"

#ifdef ALGLIB_GPU_ENABLE_CUDA

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/scan.h>
#include <thrust/device_ptr.h>

namespace alglib_gpu
{
namespace
{
inline JobStatus translate_cuda(cudaError_t err)
  {
   if(err == cudaSuccess)
      return JobStatus::Ok;
   return JobStatus::ErrorBackendFailure;
  }

struct SquareOp
  {
   __host__ __device__ double operator()(double v) const { return v * v; }
  };

class DeviceBuffer
  {
  public:
   DeviceBuffer() = default;
   DeviceBuffer(size_t bytes)
     {
      if(bytes == 0)
        {
         status_ = JobStatus::Ok;
         return;
        }
      status_ = translate_cuda(cudaMalloc(&ptr_, bytes));
     }

   ~DeviceBuffer()
     {
      if(ptr_)
         cudaFree(ptr_);
     }

   DeviceBuffer(const DeviceBuffer&) = delete;
   DeviceBuffer& operator=(const DeviceBuffer&) = delete;

   void* get() const { return ptr_; }
   JobStatus status() const { return status_; }

  private:
   void* ptr_ = nullptr;
   JobStatus status_ = JobStatus::Ok;
  };

constexpr int kBlockSize = 256;

__global__ void mask_complex_kernel(const double* spectrum,
                                    const double* mask,
                                    double*       out,
                                    std::size_t   bin_count,
                                    int           mask_is_complex,
                                    int           mode)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= bin_count)
      return;

   const double sr = spectrum[2 * idx];
   const double si = spectrum[2 * idx + 1];
   double out_r = sr;
   double out_i = si;

   if(mask_is_complex)
     {
      const double mr = mask[2 * idx];
      const double mi = mask[2 * idx + 1];
      out_r = sr * mr - si * mi;
      out_i = sr * mi + si * mr;
      if(mode == 1)
        {
         const double magnitude = sqrt(mr * mr + mi * mi);
         const double gain = (magnitude < 1.0) ? magnitude : 1.0;
         out_r *= gain;
         out_i *= gain;
        }
     }
   else
     {
      const double gain = mask[idx];
      out_r = sr * gain;
      out_i = si * gain;
     }

   out[2 * idx]     = out_r;
   out[2 * idx + 1] = out_i;
  }

__global__ void denoise_kernel(const double* spectrum,
                               double*       out,
                               std::size_t   bin_count,
                               int           method,
                               double        threshold,
                               double        beta,
                               int           iterations)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= bin_count)
      return;

   const double real = spectrum[2 * idx];
   const double imag = spectrum[2 * idx + 1];
   const double magnitude = hypot(real, imag);

   double scale = 1.0;
   if(method == 0) // hard
     {
      if(magnitude < threshold)
         scale = 0.0;
     }
   else if(method == 1) // soft
     {
      if(magnitude <= threshold)
         scale = 0.0;
      else
         scale = (magnitude - threshold) / magnitude;
     }
   else // Wiener-like
     {
      const double th2 = threshold * threshold;
      scale = magnitude / (magnitude + th2 + beta);
     }

   for(int iter = 1; iter < iterations; ++iter)
     {
      if(method == 2)
        {
         const double th2 = threshold * threshold;
         scale = scale * magnitude / (magnitude + th2 + beta);
        }
     }

   out[2 * idx]     = real * scale;
   out[2 * idx + 1] = imag * scale;
  }

__global__ void upscale_kernel(const double* spectrum,
                               double*       out,
                               std::size_t   src_bins,
                               std::size_t   dst_bins,
                               double        factor,
                               int           normalize)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= dst_bins)
      return;

   const double src_pos = static_cast<double>(idx) / factor;
   const std::size_t base_idx = static_cast<std::size_t>(floor(src_pos));
   const double frac = src_pos - static_cast<double>(base_idx);
   const std::size_t idx1 = (base_idx >= src_bins) ? src_bins - 1 : base_idx;
   const std::size_t idx2 = (base_idx + 1 >= src_bins) ? src_bins - 1 : (base_idx + 1);

   double r1 = spectrum[2 * idx1];
   double i1 = spectrum[2 * idx1 + 1];
   const double r2 = spectrum[2 * idx2];
   const double i2 = spectrum[2 * idx2 + 1];

   double out_r = r1 + (r2 - r1) * frac;
   double out_i = i1 + (i2 - i1) * frac;
   if(normalize)
     {
      out_r /= factor;
      out_i /= factor;
     }

   out[2 * idx]     = out_r;
   out[2 * idx + 1] = out_i;
  }

__global__ void downscale_kernel(const double* spectrum,
                                 double*       out,
                                 std::size_t   src_bins,
                                 std::size_t   dst_bins,
                                 double        factor,
                                 int           anti_alias)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= dst_bins)
      return;

   const double start_f = static_cast<double>(idx) * factor;
   const double end_f   = static_cast<double>(idx + 1) * factor;
   std::size_t start = static_cast<std::size_t>(floor(start_f));
   std::size_t end   = static_cast<std::size_t>(floor(end_f));
   if(end <= start)
      end = start + 1;
   if(start >= src_bins)
      start = src_bins - 1;
   if(end > src_bins)
      end = src_bins;

   double acc_r = 0.0;
   double acc_i = 0.0;
   std::size_t count = 0;
   for(std::size_t k = start; k < end; ++k)
     {
      acc_r += spectrum[2 * k];
      acc_i += spectrum[2 * k + 1];
      ++count;
     }
   if(count == 0)
      count = 1;
   double out_r = acc_r / static_cast<double>(count);
   double out_i = acc_i / static_cast<double>(count);
   if(anti_alias)
     {
      out_r *= 1.0 / factor;
      out_i *= 1.0 / factor;
     }

   out[2 * idx]     = out_r;
   out[2 * idx + 1] = out_i;
  }

__global__ void convolution_kernel(const double* lhs,
                                   const double* rhs,
                                   double*       out,
                                   std::size_t   bin_count,
                                   int           normalize)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= bin_count)
      return;

   const double xr = lhs[2 * idx];
   const double xi = lhs[2 * idx + 1];
   const double yr = rhs[2 * idx];
   const double yi = rhs[2 * idx + 1];
   double out_r = xr * yr - xi * yi;
   double out_i = xr * yi + xi * yr;
   if(normalize)
     {
      const double norm = static_cast<double>(bin_count);
      out_r /= norm;
      out_i /= norm;
     }
   out[2 * idx]     = out_r;
   out[2 * idx + 1] = out_i;
  }

__global__ void hartley_from_complex_kernel(const double* spectrum,
                                                  double*       out,
                                                  std::size_t   bins,
                                                  double        scale)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= bins)
      return;
   const std::size_t base = idx * 2;
   const double real = spectrum[base];
   const double imag = spectrum[base + 1];
   out[idx] = (real - imag) * scale;
  }

__global__ void correlation_kernel(const double* lhs,
                                   const double* rhs,
                                   double*       out,
                                   std::size_t   bin_count)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= bin_count)
      return;

   const double xr = lhs[2 * idx];
   const double xi = lhs[2 * idx + 1];
   const double yr = rhs[2 * idx];
   const double yi = rhs[2 * idx + 1];
   const double out_r = xr * yr + xi * yi;
   const double out_i = xi * yr - xr * yi;
  out[2 * idx]     = out_r;
  out[2 * idx + 1] = out_i;
 }

__global__ void resample_real_kernel(const double* input,
                                     double*       out,
                                     std::size_t   src_len,
                                     std::size_t   dst_len,
                                     double        factor)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= dst_len)
      return;

   const double src_pos = static_cast<double>(idx) / factor;
   const std::size_t base_idx = static_cast<std::size_t>(floor(src_pos));
   const double frac = src_pos - static_cast<double>(base_idx);
   const std::size_t idx1 = (base_idx >= src_len) ? src_len - 1 : base_idx;
   const std::size_t idx2 = (base_idx + 1 >= src_len) ? src_len - 1 : (base_idx + 1);
   const double v1 = input[idx1];
   const double v2 = input[idx2];
   out[idx] = v1 + (v2 - v1) * frac;
  }

__global__ void zero_pad_real_kernel(const double* input,
                                     double*       out,
                                     std::size_t   src_len,
                                     std::size_t   pad_left,
                                     std::size_t   pad_right)
  {
   const std::size_t dst_len = pad_left + src_len + pad_right;
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= dst_len)
      return;

   if(idx < pad_left || idx >= pad_left + src_len)
     {
      out[idx] = 0.0;
      return;
     }

   const std::size_t src_idx = idx - pad_left;
   out[idx] = input[src_idx];
  }

__global__ void subtract_mean_kernel(const double* input,
                                     double*       out,
                                     std::size_t   len,
                                     double        mean)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= len)
      return;
   out[idx] = input[idx] - mean;
  }

constexpr double kPi    = 3.14159265358979323846;
constexpr double kTwoPi = 6.28318530717958647692;

__global__ void compute_phase_kernel(const double2* spectrum,
                                     double*        phases,
                                     std::size_t    bins)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= bins)
      return;
   const double2 val = spectrum[idx];
   phases[idx] = atan2(val.y, val.x);
  }

__global__ void wrap_delta_kernel(const double* phases,
                                  int*          deltas,
                                  std::size_t   bins)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= bins)
      return;
   if(idx == 0)
     {
      deltas[idx] = 0;
      return;
     }
   const double delta = phases[idx] - phases[idx - 1];
   if(delta > kPi)
      deltas[idx] = -1;
   else if(delta < -kPi)
      deltas[idx] = 1;
   else
      deltas[idx] = 0;
  }

__global__ void finalize_unwrap_kernel(const double* phases,
                                       const int*    wrap_count,
                                       double*       out,
                                       std::size_t   bins)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= bins)
      return;
   out[idx] = phases[idx] + static_cast<double>(wrap_count[idx]) * kTwoPi;
  }

__global__ void instant_metrics_kernel(const double* series,
                                       const double* prefix,
                                       const double* prefix_sq,
                                       wave_pipe::InstantMetricsPayload* out,
                                       std::size_t len,
                                       int smooth_window,
                                       double epsilon)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= len)
      return;

   const int window = (smooth_window < 1) ? 1 : smooth_window;
   const int start = (idx >= static_cast<std::size_t>(window)) ? static_cast<int>(idx) - window : 0;
   int end = static_cast<int>(idx + window + 1);
   if(end > static_cast<int>(len))
      end = static_cast<int>(len);
   const int count = end - start;

   const double sum    = prefix[end - 1] - (start > 0 ? prefix[start - 1] : 0.0);
   const double sum_sq = prefix_sq[end - 1] - (start > 0 ? prefix_sq[start - 1] : 0.0);
   double mean = sum / static_cast<double>(count);
   double variance = (sum_sq / static_cast<double>(count)) - mean * mean;
   if(variance < epsilon)
      variance = epsilon;
   const double amplitude = sqrt(variance);

   wave_pipe::InstantMetricsPayload payload;
   payload.amplitude = amplitude;
   payload.energy    = amplitude * amplitude;
   payload.phase     = mean;
   payload.frequency = (idx > 0) ? (series[idx] - series[idx - 1]) : 0.0;
   payload.envelope_up   = mean + amplitude;
   payload.envelope_down = mean - amplitude;
   out[idx] = payload;
  }

__global__ void detect_transitions_kernel(const double* metrics,
                                          wave_pipe::TransitionPayload* transitions,
                                          int*       transition_count,
                                          std::size_t len,
                                          double energy_threshold,
                                          double phase_jump_threshold,
                                          int min_duration,
                                          int type_mask)
  {
   const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
   if(idx >= len)
      return;

   const double amp = fabs(metrics[idx]);
   const double prev_amp = (idx > 0) ? fabs(metrics[idx - 1]) : 0.0;
   if(amp < energy_threshold || (idx > 0 && prev_amp >= energy_threshold))
      return;

   std::size_t start = idx;
   std::size_t end   = idx;
   double last_phase = atan2(metrics[idx], (idx > 0 ? metrics[idx - 1] : metrics[idx]));

   for(std::size_t j = idx + 1; j < len; ++j)
     {
      const double current_amp = fabs(metrics[j]);
      const double phase = atan2(metrics[j], metrics[j - 1]);
      const double phase_jump = fabs(phase - last_phase);
      last_phase = phase;
      const bool duration_ok = (j > start) && (static_cast<int>(j - start) >= min_duration);
      const bool phase_trigger = phase_jump > phase_jump_threshold;
      const bool exit_amp = current_amp < energy_threshold * 0.5;
      if(duration_ok && (phase_trigger || exit_amp || j + 1 == len))
        {
         end = j;
         wave_pipe::TransitionPayload payload{};
         payload.start_index   = static_cast<int>(start);
         payload.end_index     = static_cast<int>(end);
         payload.energy_ratio  = (energy_threshold > 0.0) ? current_amp / energy_threshold : 0.0;
         payload.phase_shift   = phase_jump;
         payload.type          = type_mask;
         payload.reserved      = 0;
         const int slot = atomicAdd(transition_count, 1);
         transitions[slot] = payload;
         break;
        }
     }
  }

__global__ void remove_dc_mode1_kernel(const double* input,
                                       double*       out,
                                       std::size_t   length,
                                       double        alpha)
  {
   if(blockIdx.x != 0 || threadIdx.x != 0)
      return;
   if(length == 0)
      return;
   out[0] = 0.0;
   double prev_input = input[0];
   for(std::size_t i = 1; i < length; ++i)
     {
      const double current = input[i];
      out[i] = alpha * (out[i - 1] + current - prev_input);
      prev_input = current;
     }
  }

class CudaKernelExecutor final : public IKernelExecutor
  {
  public:
   explicit CudaKernelExecutor(const GpuConfig& cfg)
     {
      const int device = std::max(0, cfg.device_index);
      cudaSetDevice(device);
     }

   JobStatus apply_mask(const double* spectrum,
                        std::size_t   bin_count,
                        const double* mask,
                        std::size_t   mask_len,
                        bool          mask_is_complex,
                        int           mode,
                        std::vector<double>& out) override
     {
      const std::size_t spectrum_len = bin_count * 2;
      DeviceBuffer d_input(spectrum_len * sizeof(double));
      DeviceBuffer d_mask(mask_len * sizeof(double));
      DeviceBuffer d_out(spectrum_len * sizeof(double));
      if(d_input.status() != JobStatus::Ok || d_mask.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_input.get(), spectrum, spectrum_len * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaMemcpy(d_mask.get(), mask, mask_len * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((bin_count + kBlockSize - 1) / kBlockSize);
      mask_complex_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_input.get()),
                                                  static_cast<const double*>(d_mask.get()),
                                                  static_cast<double*>(d_out.get()),
                                                  bin_count,
                                                  mask_is_complex ? 1 : 0,
                                                  mode);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(spectrum_len);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), spectrum_len * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

   JobStatus denoise(const double* spectrum,
                     std::size_t   bin_count,
                     int           method,
                     double        threshold,
                     double        beta,
                     int           iterations,
                     std::vector<double>& out) override
     {
      const std::size_t spectrum_len = bin_count * 2;
      DeviceBuffer d_input(spectrum_len * sizeof(double));
      DeviceBuffer d_out(spectrum_len * sizeof(double));
      if(d_input.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_input.get(), spectrum, spectrum_len * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((bin_count + kBlockSize - 1) / kBlockSize);
      denoise_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_input.get()),
                                             static_cast<double*>(d_out.get()),
                                             bin_count,
                                             method,
                                             threshold,
                                             beta,
                                             iterations);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(spectrum_len);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), spectrum_len * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

   JobStatus upscale(const double* spectrum,
                     std::size_t   bin_count,
                     double        factor,
                     int           mode,
                     int           normalize,
                     std::vector<double>& out) override
     {
      const std::size_t dst_bins = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(static_cast<double>(bin_count) * factor)));
      const std::size_t dst_len = dst_bins * 2;
      DeviceBuffer d_input(bin_count * 2 * sizeof(double));
      DeviceBuffer d_out(dst_len * sizeof(double));
      if(d_input.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_input.get(), spectrum, bin_count * 2 * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((dst_bins + kBlockSize - 1) / kBlockSize);
      upscale_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_input.get()),
                                             static_cast<double*>(d_out.get()),
                                             bin_count,
                                             dst_bins,
                                             factor,
                                             normalize);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(dst_len);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), dst_len * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

   JobStatus downscale(const double* spectrum,
                       std::size_t   bin_count,
                       double        factor,
                       int           mode,
                       int           anti_alias,
                       std::vector<double>& out) override
     {
      (void)mode;
      const std::size_t dst_bins = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(static_cast<double>(bin_count) / factor)));
      const std::size_t dst_len = dst_bins * 2;
      DeviceBuffer d_input(bin_count * 2 * sizeof(double));
      DeviceBuffer d_out(dst_len * sizeof(double));
      if(d_input.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_input.get(), spectrum, bin_count * 2 * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((dst_bins + kBlockSize - 1) / kBlockSize);
      downscale_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_input.get()),
                                               static_cast<double*>(d_out.get()),
                                               bin_count,
                                               dst_bins,
                                               factor,
                                               anti_alias);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(dst_len);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), dst_len * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

   JobStatus convolution(const double* lhs,
                         const double* rhs,
                         std::size_t   bin_count,
                         int           normalize,
                         std::vector<double>& out) override
     {
      const std::size_t spectrum_len = bin_count * 2;
      DeviceBuffer d_lhs(spectrum_len * sizeof(double));
      DeviceBuffer d_rhs(spectrum_len * sizeof(double));
      DeviceBuffer d_out(spectrum_len * sizeof(double));
      if(d_lhs.status() != JobStatus::Ok || d_rhs.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_lhs.get(), lhs, spectrum_len * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaMemcpy(d_rhs.get(), rhs, spectrum_len * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((bin_count + kBlockSize - 1) / kBlockSize);
      convolution_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_lhs.get()),
                                                 static_cast<const double*>(d_rhs.get()),
                                                 static_cast<double*>(d_out.get()),
                                                 bin_count,
                                                 normalize);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(spectrum_len);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), spectrum_len * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

   JobStatus correlation(const double* lhs,
                         const double* rhs,
                         std::size_t   bin_count,
                         std::vector<double>& out) override
     {
      const std::size_t spectrum_len = bin_count * 2;
      DeviceBuffer d_lhs(spectrum_len * sizeof(double));
      DeviceBuffer d_rhs(spectrum_len * sizeof(double));
      DeviceBuffer d_out(spectrum_len * sizeof(double));
      if(d_lhs.status() != JobStatus::Ok || d_rhs.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_lhs.get(), lhs, spectrum_len * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaMemcpy(d_rhs.get(), rhs, spectrum_len * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((bin_count + kBlockSize - 1) / kBlockSize);
      correlation_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_lhs.get()),
                                                 static_cast<const double*>(d_rhs.get()),
                                                 static_cast<double*>(d_out.get()),
                                                 bin_count);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(spectrum_len);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), spectrum_len * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

   JobStatus resample(const double* input,
                      std::size_t   length,
                      double        factor,
                      double        cutoff,
                      int           method,
                      std::vector<double>& out) override
     {
      (void)cutoff;
      (void)method;
      if(length == 0)
        {
         out.clear();
         return JobStatus::Ok;
        }
      if(factor <= 1.0000001)
        {
         out.assign(input, input + length);
         return JobStatus::Ok;
        }

      const std::size_t dst_len = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(static_cast<double>(length) * factor)));
      DeviceBuffer d_input(length * sizeof(double));
      DeviceBuffer d_out(dst_len * sizeof(double));
      if(d_input.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_input.get(), input, length * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((dst_len + kBlockSize - 1) / kBlockSize);
      resample_real_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_input.get()),
                                                   static_cast<double*>(d_out.get()),
                                                   length,
                                                   dst_len,
                                                   factor);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(dst_len);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), dst_len * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

   JobStatus zero_pad(const double* input,
                      std::size_t   length,
                      std::size_t   pad_left,
                      std::size_t   pad_right,
                      std::vector<double>& out) override
     {
      const std::size_t dst_len = pad_left + length + pad_right;
      DeviceBuffer d_input(length * sizeof(double));
      DeviceBuffer d_out(dst_len * sizeof(double));
      if(d_input.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_input.get(), input, length * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((dst_len + kBlockSize - 1) / kBlockSize);
      zero_pad_real_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_input.get()),
                                                   static_cast<double*>(d_out.get()),
                                                   length,
                                                   pad_left,
                                                   pad_right);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(dst_len);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), dst_len * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

  JobStatus remove_dc(const double* input,
                       std::size_t   length,
                       int           mode,
                       double        alpha,
                       std::vector<double>& out) override
     {
      if(length == 0)
        {
         out.clear();
         return JobStatus::Ok;
        }

      DeviceBuffer d_input(length * sizeof(double));
      DeviceBuffer d_out(length * sizeof(double));
      if(d_input.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_input.get(), input, length * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      if(mode == 0)
        {
         double sum = 0.0;
         for(std::size_t i = 0; i < length; ++i)
            sum += input[i];
         const double mean = sum / static_cast<double>(length);

         const int blocks = static_cast<int>((length + kBlockSize - 1) / kBlockSize);
         subtract_mean_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_input.get()),
                                                      static_cast<double*>(d_out.get()),
                                                      length,
                                                      mean);
         status = translate_cuda(cudaGetLastError());
         if(status != JobStatus::Ok)
            return status;
         status = translate_cuda(cudaDeviceSynchronize());
         if(status != JobStatus::Ok)
            return status;
        }
      else if(mode == 1)
        {
         remove_dc_mode1_kernel<<<1, 1>>>(static_cast<const double*>(d_input.get()),
                                           static_cast<double*>(d_out.get()),
                                           length,
                                           alpha);
         status = translate_cuda(cudaGetLastError());
         if(status != JobStatus::Ok)
            return status;
         status = translate_cuda(cudaDeviceSynchronize());
         if(status != JobStatus::Ok)
            return status;
        }
      else
        {
         return JobStatus::ErrorInvalidArgs;
        }

      out.resize(length);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), length * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

  JobStatus phase_unwrap(const double* spectrum,
                         std::size_t   bin_count,
                         int           method,
                         std::vector<double>& out) override
     {
      (void)method;
      if(bin_count == 0)
        {
         out.clear();
         return JobStatus::Ok;
        }

      thrust::device_vector<double2> d_spectrum(bin_count);
      auto status = translate_cuda(cudaMemcpy(thrust::raw_pointer_cast(d_spectrum.data()),
                                              spectrum,
                                              bin_count * sizeof(double2),
                                              cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      thrust::device_vector<double> d_phases(bin_count);
      const int blocks = static_cast<int>((bin_count + kBlockSize - 1) / kBlockSize);
      compute_phase_kernel<<<blocks, kBlockSize>>>(thrust::raw_pointer_cast(d_spectrum.data()),
                                                  thrust::raw_pointer_cast(d_phases.data()),
                                                  bin_count);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      thrust::device_vector<int> d_deltas(bin_count);
      wrap_delta_kernel<<<blocks, kBlockSize>>>(thrust::raw_pointer_cast(d_phases.data()),
                                               thrust::raw_pointer_cast(d_deltas.data()),
                                               bin_count);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      thrust::device_vector<int> d_wrap(bin_count);
      thrust::exclusive_scan(d_deltas.begin(), d_deltas.end(), d_wrap.begin(), 0);

      thrust::device_vector<double> d_unwrapped(bin_count);
      finalize_unwrap_kernel<<<blocks, kBlockSize>>>(thrust::raw_pointer_cast(d_phases.data()),
                                                     thrust::raw_pointer_cast(d_wrap.data()),
                                                     thrust::raw_pointer_cast(d_unwrapped.data()),
                                                     bin_count);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(bin_count);
      status = translate_cuda(cudaMemcpy(out.data(), thrust::raw_pointer_cast(d_unwrapped.data()),
                                         bin_count * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

  JobStatus instant_metrics(const double* time_series,
                             std::size_t   length,
                             int           smooth_window,
                             double        epsilon,
                             std::vector<wave_pipe::InstantMetricsPayload>& instant_out) override
     {
      if(length == 0)
        {
         instant_out.clear();
         return JobStatus::Ok;
        }

      thrust::device_vector<double> d_series(length);
      auto status = translate_cuda(cudaMemcpy(thrust::raw_pointer_cast(d_series.data()),
                                              time_series,
                                              length * sizeof(double),
                                              cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      thrust::device_vector<double> d_prefix(length);
      thrust::inclusive_scan(d_series.begin(), d_series.end(), d_prefix.begin());

      thrust::device_vector<double> d_prefix_sq(length);
      thrust::transform(d_series.begin(), d_series.end(), d_prefix_sq.begin(), SquareOp{});
      thrust::inclusive_scan(d_prefix_sq.begin(), d_prefix_sq.end(), d_prefix_sq.begin());

      thrust::device_vector<wave_pipe::InstantMetricsPayload> d_metrics(length);
      const int blocks = static_cast<int>((length + kBlockSize - 1) / kBlockSize);
      instant_metrics_kernel<<<blocks, kBlockSize>>>(thrust::raw_pointer_cast(d_series.data()),
                                                     thrust::raw_pointer_cast(d_prefix.data()),
                                                     thrust::raw_pointer_cast(d_prefix_sq.data()),
                                                     thrust::raw_pointer_cast(d_metrics.data()),
                                                     length,
                                                     smooth_window,
                                                     epsilon);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      instant_out.resize(length);
      status = translate_cuda(cudaMemcpy(instant_out.data(), thrust::raw_pointer_cast(d_metrics.data()),
                                         length * sizeof(wave_pipe::InstantMetricsPayload), cudaMemcpyDeviceToHost));
      return status;
     }

  JobStatus detect_transitions(const double* metrics,
                               std::size_t   length,
                               double        energy_threshold,
                               double        phase_jump_threshold,
                               int           min_duration,
                               int           type_mask,
                               std::vector<wave_pipe::TransitionPayload>& transitions) override
     {
      if(length == 0)
        {
         transitions.clear();
         return JobStatus::Ok;
        }

      thrust::device_vector<double> d_metrics(length);
      auto status = translate_cuda(cudaMemcpy(thrust::raw_pointer_cast(d_metrics.data()),
                                              metrics,
                                              length * sizeof(double),
                                              cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      thrust::device_vector<wave_pipe::TransitionPayload> d_transitions(length);
      thrust::device_vector<int> d_count(1, 0);

      const int blocks = static_cast<int>((length + kBlockSize - 1) / kBlockSize);
      detect_transitions_kernel<<<blocks, kBlockSize>>>(thrust::raw_pointer_cast(d_metrics.data()),
                                                       thrust::raw_pointer_cast(d_transitions.data()),
                                                       thrust::raw_pointer_cast(d_count.data()),
                                                       length,
                                                       energy_threshold,
                                                       phase_jump_threshold,
                                                       min_duration,
                                                       type_mask);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      int host_count = 0;
      status = translate_cuda(cudaMemcpy(&host_count, thrust::raw_pointer_cast(d_count.data()), sizeof(int), cudaMemcpyDeviceToHost));
      if(status != JobStatus::Ok)
         return status;
      host_count = std::min(host_count, static_cast<int>(length));

      transitions.resize(host_count);
      if(host_count > 0)
        status = translate_cuda(cudaMemcpy(transitions.data(), thrust::raw_pointer_cast(d_transitions.data()),
                                           host_count * sizeof(wave_pipe::TransitionPayload), cudaMemcpyDeviceToHost));
      return status;
     }

  JobStatus hartley_from_complex(const double* spectrum,
                                 std::size_t   bin_count,
                                 double        scale,
                                 std::vector<double>& out) override
     {
      if(bin_count == 0)
        {
         out.clear();
         return JobStatus::Ok;
        }

      const std::size_t spectrum_len = bin_count * 2;
      DeviceBuffer d_input(spectrum_len * sizeof(double));
      DeviceBuffer d_out(bin_count * sizeof(double));
      if(d_input.status() != JobStatus::Ok || d_out.status() != JobStatus::Ok)
         return JobStatus::ErrorBackendFailure;

      auto status = translate_cuda(cudaMemcpy(d_input.get(), spectrum, spectrum_len * sizeof(double), cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
         return status;

      const int blocks = static_cast<int>((bin_count + kBlockSize - 1) / kBlockSize);
      hartley_from_complex_kernel<<<blocks, kBlockSize>>>(static_cast<const double*>(d_input.get()),
                                                          static_cast<double*>(d_out.get()),
                                                          bin_count,
                                                          scale);
      status = translate_cuda(cudaGetLastError());
      if(status != JobStatus::Ok)
         return status;
      status = translate_cuda(cudaDeviceSynchronize());
      if(status != JobStatus::Ok)
         return status;

      out.resize(bin_count);
      status = translate_cuda(cudaMemcpy(out.data(), d_out.get(), bin_count * sizeof(double), cudaMemcpyDeviceToHost));
      return status;
     }

  };
} // namespace

std::unique_ptr<IKernelExecutor> CreateCudaKernelExecutor(const GpuConfig& cfg)
  {
   try
     {
      return std::make_unique<CudaKernelExecutor>(cfg);
     }
   catch(...)
     {
      return nullptr;
     }
  }

} // namespace alglib_gpu

#endif // ALGLIB_GPU_ENABLE_CUDA
