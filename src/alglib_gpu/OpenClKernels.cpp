#include "GpuBackend.h"

#ifdef ALGLIB_GPU_ENABLE_OPENCL

// Normalizar alvo e suprimir APIs 1.2 como deprecadas para build limpo com /WX
#ifndef CL_TARGET_OPENCL_VERSION
#define CL_TARGET_OPENCL_VERSION 120
#endif
#ifndef CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#endif
#include <CL/cl.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace alglib_gpu
{
namespace
{
inline JobStatus translate_cl(cl_int err)
  {
   return (err == CL_SUCCESS) ? JobStatus::Ok : JobStatus::ErrorBackendFailure;
  }

constexpr const char* kKernelSource = R"(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

typedef struct {
   double amplitude;
   double frequency;
   double energy;
   double phase;
   double envelope_up;
   double envelope_down;
} InstantMetricsPayload;

typedef struct {
   int    start_index;
   int    end_index;
   double energy_ratio;
   double phase_shift;
   int    type;
   int    reserved;
} TransitionPayload;

__kernel void mask_complex_kernel(__global const double2* spectrum,
                                  __global const double2* mask,
                                  __global double2*       out,
                                  const int               mode)
{
   const int idx = get_global_id(0);
   const double2 s = spectrum[idx];
   const double2 m = mask[idx];
   double2 res;
   res.x = s.x * m.x - s.y * m.y;
   res.y = s.x * m.y + s.y * m.x;
   if(mode == 1)
     {
      const double magnitude = hypot(m.x, m.y);
      const double gain = magnitude < 1.0 ? magnitude : 1.0;
      res *= gain;
     }
   out[idx] = res;
}

__kernel void mask_real_kernel(__global const double2* spectrum,
                               __global const double*  mask,
                               __global double2*       out,
                               const int               mode)
{
   const int idx = get_global_id(0);
   const double2 s = spectrum[idx];
   double gain = mask[idx];
   if(mode == 1)
      gain = gain < 1.0 ? gain : 1.0;
   double2 res = (double2)(s.x * gain, s.y * gain);
   out[idx] = res;
}

__kernel void denoise_kernel(__global const double2* spectrum,
                             __global double2*       out,
                             const int               method,
                             const double            threshold,
                             const double            beta,
                             const int               iterations)
{
   const int idx = get_global_id(0);
   const double2 val = spectrum[idx];
   const double magnitude = hypot(val.x, val.y);
   double scale = 1.0;
   if(method == 0)
     {
      if(magnitude < threshold)
         scale = 0.0;
     }
   else if(method == 1)
     {
      if(magnitude <= threshold)
         scale = 0.0;
      else
         scale = (magnitude - threshold) / magnitude;
     }
   else
     {
      const double th2 = threshold * threshold;
      scale = magnitude / (magnitude + th2 + beta);
      for(int iter = 1; iter < iterations; ++iter)
        scale = scale * magnitude / (magnitude + th2 + beta);
     }
   out[idx] = (double2)(val.x * scale, val.y * scale);
}

__kernel void upscale_kernel(__global const double2* spectrum,
                             __global double2*       out,
                             const int               src_bins,
                             const double            factor,
                             const int               normalize)
{
   const int idx = get_global_id(0);
   const double src_pos = (double)idx / factor;
   const int base_idx = (int)floor(src_pos);
   const double frac = src_pos - (double)base_idx;
   const int idx1 = base_idx >= src_bins ? src_bins - 1 : base_idx;
   const int idx2 = base_idx + 1 >= src_bins ? src_bins - 1 : base_idx + 1;
   const double2 a = spectrum[idx1];
   const double2 b = spectrum[idx2];
   double2 res = (double2)(a.x + (b.x - a.x) * frac,
                           a.y + (b.y - a.y) * frac);
   if(normalize)
     res /= factor;
   out[idx] = res;
}

__kernel void downscale_kernel(__global const double2* spectrum,
                               __global double2*       out,
                               const int               src_bins,
                               const double            factor,
                               const int               anti_alias)
{
   const int idx = get_global_id(0);
   const double start_f = (double)idx * factor;
   const double end_f   = (double)(idx + 1) * factor;
   int start = (int)floor(start_f);
   int end   = (int)floor(end_f);
   if(end <= start)
      end = start + 1;
   if(start >= src_bins)
      start = src_bins - 1;
   if(end > src_bins)
      end = src_bins;

   double2 acc = (double2)(0.0, 0.0);
   int count = 0;
   for(int k = start; k < end; ++k)
     {
      acc += spectrum[k];
      ++count;
     }
   if(count == 0)
      count = 1;
   acc /= (double)count;
   if(anti_alias)
      acc /= factor;
   out[idx] = acc;
}

__kernel void convolution_kernel(__global const double2* lhs,
                                 __global const double2* rhs,
                                 __global double2*       out,
                                 const double            norm)
{
   const int idx = get_global_id(0);
   const double2 a = lhs[idx];
   const double2 b = rhs[idx];
   double2 res = (double2)((a.x * b.x - a.y * b.y) / norm,
                           (a.x * b.y + a.y * b.x) / norm);
   out[idx] = res;
}

__kernel void correlation_kernel(__global const double2* lhs,
                                 __global const double2* rhs,
                                 __global double2*       out)
{
   const int idx = get_global_id(0);
   const double2 a = lhs[idx];
   const double2 b = rhs[idx];
   double2 res = (double2)(a.x * b.x + a.y * b.y,
                           a.y * b.x - a.x * b.y);
   out[idx] = res;
}

__kernel void resample_kernel(__global const double* input,
                              __global double*       out,
                              const int              src_len,
                              const double           factor)
{
   const int idx = get_global_id(0);
   const double src_pos = (double)idx / factor;
   const int base_idx = (int)floor(src_pos);
   const double frac = src_pos - (double)base_idx;
   const int idx1 = base_idx >= src_len ? src_len - 1 : base_idx;
   const int idx2 = base_idx + 1 >= src_len ? src_len - 1 : base_idx + 1;
   const double v1 = input[idx1];
   const double v2 = input[idx2];
   out[idx] = v1 + (v2 - v1) * frac;
}

__kernel void zero_pad_kernel(__global const double* input,
                              __global double*       out,
                              const int              src_len,
                              const int              pad_left)
{
   const int idx = get_global_id(0);
   if(idx < pad_left || idx >= pad_left + src_len)
      out[idx] = 0.0;
   else
      out[idx] = input[idx - pad_left];
}

__kernel void subtract_mean_kernel(__global const double* input,
                                   __global double*       out,
                                   const double           mean)
{
   const int idx = get_global_id(0);
   out[idx] = input[idx] - mean;
}

__kernel void remove_dc_mode1_kernel(__global const double* input,
                                     __global double*       out,
                                     const int              length,
                                     const double           alpha)
{
   if(get_global_id(0) != 0)
      return;
   if(length <= 0)
      return;
   out[0] = 0.0;
   double prev_input = input[0];
   for(int i = 1; i < length; ++i)
     {
      const double current = input[i];
      out[i] = alpha * (out[i - 1] + current - prev_input);
      prev_input = current;
     }
}

__kernel void compute_phase_kernel(__global const double2* spectrum,
                                   __global double*        phases,
                                   const int               bins)
{
   const int idx = get_global_id(0);
   if(idx >= bins)
      return;
   const double2 val = spectrum[idx];
   phases[idx] = atan2(val.y, val.x);
}

__kernel void wrap_delta_kernel(__global const double* phases,
                                __global int*          deltas,
                                const int              bins)
{
   const int idx = get_global_id(0);
   if(idx >= bins)
      return;
   if(idx == 0)
     {
      deltas[idx] = 0;
      return;
     }
   const double prev = phases[idx - 1];
   const double current = phases[idx];
   const double delta = current - prev;
   const double pi = 3.14159265358979323846;
   if(delta > pi)
      deltas[idx] = -1;
   else if(delta < -pi)
      deltas[idx] = 1;
   else
      deltas[idx] = 0;
}

__kernel void exclusive_scan_kernel(__global const int* deltas,
                                    __global int*       wrap_count,
                                    __global int*       block_totals,
                                    const int           bins,
                                    __local int*        temp)
{
   const int lid        = get_local_id(0);
   const int local_size = get_local_size(0);
   const int group      = get_group_id(0);
   const int index      = group * local_size + lid;

   int value = (index < bins) ? deltas[index] : 0;
   temp[lid] = value;
   barrier(CLK_LOCAL_MEM_FENCE);

   for(int offset = 1; offset < local_size; offset <<= 1)
     {
      int prev = (lid >= offset) ? temp[lid - offset] : 0;
      barrier(CLK_LOCAL_MEM_FENCE);
      temp[lid] += prev;
      barrier(CLK_LOCAL_MEM_FENCE);
     }

   if(index < bins)
     wrap_count[index] = temp[lid] - value;

   if(lid == local_size - 1 || index == bins - 1)
      block_totals[group] = temp[lid];
}

__kernel void exclusive_scan_add_kernel(__global int*       wrap_count,
                                        __global const int* block_offsets,
                                        const int           bins)
{
   const int idx = get_global_id(0);
   if(idx >= bins)
      return;
   const int group = get_group_id(0);
   wrap_count[idx] += block_offsets[group];
}

__kernel void finalize_unwrap_kernel(__global const double* phases,
                                     __global const int*    wrap_count,
                                     __global double*       out,
                                     const int              bins)
{
   const int idx = get_global_id(0);
   if(idx >= bins)
      return;
   const double two_pi = 6.28318530717958647692;
   out[idx] = phases[idx] + (double)wrap_count[idx] * two_pi;
}

__kernel void instant_metrics_kernel(__global const double* series,
                                     __global InstantMetricsPayload* out,
                                     const int              len,
                                     const int              smooth_window,
                                     const double           epsilon)
{
   const int idx = get_global_id(0);
   if(idx >= len)
      return;

   const int window = smooth_window < 1 ? 1 : smooth_window;
   int start = idx - window;
   if(start < 0)
      start = 0;
   int end = idx + window + 1;
   if(end > len)
      end = len;
   const int count = end - start;

   double sum = 0.0;
   double sum_sq = 0.0;
   for(int i = start; i < end; ++i)
     {
      const double v = series[i];
      sum    += v;
      sum_sq += v * v;
     }
   const double mean = sum / (double)count;
   double variance = (sum_sq / (double)count) - mean * mean;
   if(variance < epsilon)
      variance = epsilon;
   const double amplitude = sqrt(variance);
   const double current = series[idx];
   const double prev = (idx > 0) ? series[idx - 1] : current;
   const double frequency = current - prev;

   InstantMetricsPayload payload;
   payload.amplitude     = amplitude;
   payload.frequency     = frequency;
   payload.energy        = amplitude * amplitude;
   payload.phase         = mean;
   payload.envelope_up   = mean + amplitude;
   payload.envelope_down = mean - amplitude;
   out[idx] = payload;
}

__kernel void hartley_from_complex_kernel(__global const double2* spectrum,
                                                 __global double*       out,
                                                 const double           scale)
{
   const int idx = get_global_id(0);
   out[idx] = (spectrum[idx].x - spectrum[idx].y) * scale;
}

__kernel void detect_transitions_kernel(__global const double* metrics,
                                        __global TransitionPayload* transitions,
                                        volatile __global int*      transition_count,
                                        const int                   len,
                                        const double                energy_threshold,
                                        const double                phase_jump_threshold,
                                        const int                   min_duration,
                                        const int                   type_mask)
{
   const int idx = get_global_id(0);
   if(idx >= len)
      return;

   const double current = metrics[idx];
   const double amp = fabs(current);
   const double prev_amp = (idx > 0) ? fabs(metrics[idx - 1]) : 0.0;
   if(amp < energy_threshold || (idx > 0 && prev_amp >= energy_threshold))
      return;

   int start = idx;
   int end   = idx;
   double last_phase = atan2(current, (idx > 0) ? metrics[idx - 1] : current);

   for(int j = idx + 1; j < len; ++j)
     {
      const double value      = metrics[j];
      const double prev_value = metrics[j - 1];
      const double current_amp = fabs(value);
      const double phase = atan2(value, prev_value);
      const double phase_jump = fabs(phase - last_phase);
      last_phase = phase;

      if(current_amp < energy_threshold * 0.5)
         break;

      end = j;

      if((j - start) >= min_duration &&
         (phase_jump > phase_jump_threshold || j + 1 == len))
        {
         TransitionPayload payload;
         payload.start_index  = start;
         payload.end_index    = end;
         payload.energy_ratio = (energy_threshold > 0.0) ? current_amp / energy_threshold : 0.0;
         payload.phase_shift  = phase_jump;
         payload.type         = type_mask;
         payload.reserved     = 0;

         int slot = atomic_inc(transition_count);
         if(slot < len)
            transitions[slot] = payload;
         break;
        }
     }
}
)";

class OpenClKernelExecutor final : public IKernelExecutor
  {
  public:
   explicit OpenClKernelExecutor(const GpuConfig& cfg)
     {
      cl_uint platform_count = 0;
      cl_int err = clGetPlatformIDs(0, nullptr, &platform_count);
      if(err != CL_SUCCESS || platform_count == 0)
         throw std::runtime_error("Nenhuma plataforma OpenCL disponível");

      std::vector<cl_platform_id> platforms(platform_count);
      clGetPlatformIDs(platform_count, platforms.data(), nullptr);

      struct DeviceEntry { cl_platform_id platform; cl_device_id device; };
      std::vector<DeviceEntry> devices;
      for(auto platform : platforms)
        {
         cl_uint device_count = 0;
         err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &device_count);
         if(err != CL_SUCCESS || device_count == 0)
            continue;

         std::vector<cl_device_id> local(device_count);
         err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, device_count, local.data(), nullptr);
         if(err != CL_SUCCESS)
            continue;

         for(auto device : local)
            devices.push_back({platform, device});
        }

      if(devices.empty())
         throw std::runtime_error("Nenhum dispositivo OpenCL GPU encontrado");

      const size_t index = std::min(static_cast<size_t>(std::max(0, cfg.device_index)), devices.size() - 1);
      platform_ = devices[index].platform;
      device_   = devices[index].device;

      context_ = clCreateContext(nullptr, 1, &device_, nullptr, nullptr, &err);
      if(err != CL_SUCCESS)
         throw std::runtime_error("Falha ao criar contexto OpenCL");

      queue_ = clCreateCommandQueue(context_, device_, 0, &err);
      if(err != CL_SUCCESS)
         throw std::runtime_error("Falha ao criar command queue OpenCL");

      const char* src = kKernelSource;
      const size_t len = std::strlen(kKernelSource);
      program_ = clCreateProgramWithSource(context_, 1, &src, &len, &err);
      if(err != CL_SUCCESS)
         throw std::runtime_error("Falha ao criar programa OpenCL");

      err = clBuildProgram(program_, 1, &device_, nullptr, nullptr, nullptr);
      if(err != CL_SUCCESS)
        {
         size_t log_size = 0;
         clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
         std::string log(log_size, '\0');
         clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_LOG, log_size, log.data(), nullptr);
         throw std::runtime_error("Erro ao compilar kernels OpenCL:\n" + log);
        }

      mask_complex_kernel_ = clCreateKernel(program_, "mask_complex_kernel", &err);
      mask_real_kernel_    = clCreateKernel(program_, "mask_real_kernel", &err);
      denoise_kernel_      = clCreateKernel(program_, "denoise_kernel", &err);
      upscale_kernel_      = clCreateKernel(program_, "upscale_kernel", &err);
      downscale_kernel_    = clCreateKernel(program_, "downscale_kernel", &err);
      convolution_kernel_  = clCreateKernel(program_, "convolution_kernel", &err);
      correlation_kernel_  = clCreateKernel(program_, "correlation_kernel", &err);
      resample_kernel_     = clCreateKernel(program_, "resample_kernel", &err);
      zero_pad_kernel_     = clCreateKernel(program_, "zero_pad_kernel", &err);
      subtract_mean_kernel_= clCreateKernel(program_, "subtract_mean_kernel", &err);
      hartley_kernel_      = clCreateKernel(program_, "hartley_from_complex_kernel", &err);
      remove_dc_mode1_kernel_ = clCreateKernel(program_, "remove_dc_mode1_kernel", &err);
      compute_phase_kernel_   = clCreateKernel(program_, "compute_phase_kernel", &err);
      wrap_delta_kernel_      = clCreateKernel(program_, "wrap_delta_kernel", &err);
      exclusive_scan_kernel_  = clCreateKernel(program_, "exclusive_scan_kernel", &err);
      exclusive_scan_add_kernel_ = clCreateKernel(program_, "exclusive_scan_add_kernel", &err);
      finalize_unwrap_kernel_ = clCreateKernel(program_, "finalize_unwrap_kernel", &err);
      instant_metrics_kernel_ = clCreateKernel(program_, "instant_metrics_kernel", &err);
      detect_transitions_kernel_ = clCreateKernel(program_, "detect_transitions_kernel", &err);
     }

   ~OpenClKernelExecutor() override
     {
      if(mask_complex_kernel_) clReleaseKernel(mask_complex_kernel_);
      if(mask_real_kernel_)    clReleaseKernel(mask_real_kernel_);
      if(denoise_kernel_)      clReleaseKernel(denoise_kernel_);
      if(upscale_kernel_)      clReleaseKernel(upscale_kernel_);
      if(downscale_kernel_)    clReleaseKernel(downscale_kernel_);
      if(convolution_kernel_)  clReleaseKernel(convolution_kernel_);
      if(correlation_kernel_)  clReleaseKernel(correlation_kernel_);
      if(resample_kernel_)     clReleaseKernel(resample_kernel_);
      if(zero_pad_kernel_)     clReleaseKernel(zero_pad_kernel_);
      if(subtract_mean_kernel_)clReleaseKernel(subtract_mean_kernel_);
      if(remove_dc_mode1_kernel_) clReleaseKernel(remove_dc_mode1_kernel_);
      if(compute_phase_kernel_)   clReleaseKernel(compute_phase_kernel_);
      if(wrap_delta_kernel_)      clReleaseKernel(wrap_delta_kernel_);
      if(exclusive_scan_kernel_)  clReleaseKernel(exclusive_scan_kernel_);
      if(exclusive_scan_add_kernel_) clReleaseKernel(exclusive_scan_add_kernel_);
      if(finalize_unwrap_kernel_) clReleaseKernel(finalize_unwrap_kernel_);
      if(instant_metrics_kernel_) clReleaseKernel(instant_metrics_kernel_);
      if(detect_transitions_kernel_) clReleaseKernel(detect_transitions_kernel_);
      if(hartley_kernel_)      clReleaseKernel(hartley_kernel_);
      if(program_)             clReleaseProgram(program_);
      if(queue_)               clReleaseCommandQueue(queue_);
      if(context_)             clReleaseContext(context_);
     }

   JobStatus apply_mask(const double* spectrum,
                        std::size_t   bin_count,
                        const double* mask,
                        std::size_t   mask_len,
                        bool          mask_is_complex,
                        int           mode,
                        std::vector<double>& out) override
     {
      const size_t spectrum_bytes = bin_count * sizeof(cl_double2);
      std::vector<cl_double2> spectrum_vec(bin_count);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         spectrum_vec[i].s[0] = spectrum[2 * i];
         spectrum_vec[i].s[1] = spectrum[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem spectrum_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              spectrum_bytes, spectrum_vec.data(), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      cl_mem mask_buffer = nullptr;
      if(mask_is_complex)
        {
         std::vector<cl_double2> mask_vec(bin_count);
         for(std::size_t i = 0; i < bin_count; ++i)
           {
            mask_vec[i].s[0] = mask[2 * i];
            mask_vec[i].s[1] = mask[2 * i + 1];
           }
         mask_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                      spectrum_bytes, mask_vec.data(), &err);
        }
      else
        {
         mask_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                      mask_len * sizeof(double), const_cast<double*>(mask), &err);
        }
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         if(mask_buffer) clReleaseMemObject(mask_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, spectrum_bytes, nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(mask_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      cl_kernel kernel = mask_is_complex ? mask_complex_kernel_ : mask_real_kernel_;
      err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &spectrum_buffer);
      err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &mask_buffer);
      err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &output_buffer);
      err |= clSetKernelArg(kernel, 3, sizeof(int), &mode);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(mask_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = bin_count;
      err = clEnqueueNDRangeKernel(queue_, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(mask_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<cl_double2> output_vec(bin_count);
      err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, spectrum_bytes,
                                output_vec.data(), 0, nullptr, nullptr);

      clReleaseMemObject(spectrum_buffer);
      clReleaseMemObject(mask_buffer);
      clReleaseMemObject(output_buffer);

      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      out.resize(bin_count * 2);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         out[2 * i]     = output_vec[i].s[0];
         out[2 * i + 1] = output_vec[i].s[1];
        }
      return JobStatus::Ok;
     }

   JobStatus denoise(const double* spectrum,
                     std::size_t   bin_count,
                     int           method,
                     double        threshold,
                     double        beta,
                     int           iterations,
                     std::vector<double>& out) override
     {
      const size_t bytes = bin_count * sizeof(cl_double2);
      std::vector<cl_double2> spectrum_vec(bin_count);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         spectrum_vec[i].s[0] = spectrum[2 * i];
         spectrum_vec[i].s[1] = spectrum[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem spectrum_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              bytes, spectrum_vec.data(), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, bytes, nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(denoise_kernel_, 0, sizeof(cl_mem), &spectrum_buffer);
      err |= clSetKernelArg(denoise_kernel_, 1, sizeof(cl_mem), &output_buffer);
      err |= clSetKernelArg(denoise_kernel_, 2, sizeof(int), &method);
      err |= clSetKernelArg(denoise_kernel_, 3, sizeof(double), &threshold);
      err |= clSetKernelArg(denoise_kernel_, 4, sizeof(double), &beta);
      err |= clSetKernelArg(denoise_kernel_, 5, sizeof(int), &iterations);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = bin_count;
      err = clEnqueueNDRangeKernel(queue_, denoise_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<cl_double2> output_vec(bin_count);
      err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, bytes, output_vec.data(), 0, nullptr, nullptr);

      clReleaseMemObject(spectrum_buffer);
      clReleaseMemObject(output_buffer);

      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      out.resize(bin_count * 2);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         out[2 * i]     = output_vec[i].s[0];
         out[2 * i + 1] = output_vec[i].s[1];
        }
      return JobStatus::Ok;
     }

   JobStatus upscale(const double* spectrum,
                     std::size_t   bin_count,
                     double        factor,
                     int           mode,
                     int           normalize,
                     std::vector<double>& out) override
     {
      (void)mode;
      const size_t dst_bins = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(static_cast<double>(bin_count) * factor)));
      const size_t src_bytes = bin_count * sizeof(cl_double2);
      std::vector<cl_double2> spectrum_vec(bin_count);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         spectrum_vec[i].s[0] = spectrum[2 * i];
         spectrum_vec[i].s[1] = spectrum[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem spectrum_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              src_bytes, spectrum_vec.data(), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, dst_bins * sizeof(cl_double2), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const int src_bins_int = static_cast<int>(bin_count);
      err  = clSetKernelArg(upscale_kernel_, 0, sizeof(cl_mem), &spectrum_buffer);
      err |= clSetKernelArg(upscale_kernel_, 1, sizeof(cl_mem), &output_buffer);
      err |= clSetKernelArg(upscale_kernel_, 2, sizeof(int), &src_bins_int);
      err |= clSetKernelArg(upscale_kernel_, 3, sizeof(double), &factor);
      err |= clSetKernelArg(upscale_kernel_, 4, sizeof(int), &normalize);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = dst_bins;
      err = clEnqueueNDRangeKernel(queue_, upscale_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<cl_double2> output_vec(dst_bins);
      err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, dst_bins * sizeof(cl_double2), output_vec.data(), 0, nullptr, nullptr);

      clReleaseMemObject(spectrum_buffer);
      clReleaseMemObject(output_buffer);

      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      out.resize(dst_bins * 2);
      for(std::size_t i = 0; i < dst_bins; ++i)
        {
         out[2 * i]     = output_vec[i].s[0];
         out[2 * i + 1] = output_vec[i].s[1];
        }
      return JobStatus::Ok;
     }

   JobStatus downscale(const double* spectrum,
                       std::size_t   bin_count,
                       double        factor,
                       int           mode,
                       int           anti_alias,
                       std::vector<double>& out) override
     {
      (void)mode;
      const size_t dst_bins = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(static_cast<double>(bin_count) / factor)));
      const size_t src_bytes = bin_count * sizeof(cl_double2);
      std::vector<cl_double2> spectrum_vec(bin_count);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         spectrum_vec[i].s[0] = spectrum[2 * i];
         spectrum_vec[i].s[1] = spectrum[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem spectrum_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              src_bytes, spectrum_vec.data(), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, dst_bins * sizeof(cl_double2), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const int src_bins_int = static_cast<int>(bin_count);
      err  = clSetKernelArg(downscale_kernel_, 0, sizeof(cl_mem), &spectrum_buffer);
      err |= clSetKernelArg(downscale_kernel_, 1, sizeof(cl_mem), &output_buffer);
      err |= clSetKernelArg(downscale_kernel_, 2, sizeof(int), &src_bins_int);
      err |= clSetKernelArg(downscale_kernel_, 3, sizeof(double), &factor);
      err |= clSetKernelArg(downscale_kernel_, 4, sizeof(int), &anti_alias);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = dst_bins;
      err = clEnqueueNDRangeKernel(queue_, downscale_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<cl_double2> output_vec(dst_bins);
      err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, dst_bins * sizeof(cl_double2), output_vec.data(), 0, nullptr, nullptr);

      clReleaseMemObject(spectrum_buffer);
      clReleaseMemObject(output_buffer);

      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      out.resize(dst_bins * 2);
      for(std::size_t i = 0; i < dst_bins; ++i)
        {
         out[2 * i]     = output_vec[i].s[0];
         out[2 * i + 1] = output_vec[i].s[1];
        }
      return JobStatus::Ok;
     }

   JobStatus convolution(const double* lhs,
                         const double* rhs,
                         std::size_t   bin_count,
                         int           normalize,
                         std::vector<double>& out) override
     {
      const double norm = normalize ? static_cast<double>(bin_count) : 1.0;
      const size_t bytes = bin_count * sizeof(cl_double2);
      std::vector<cl_double2> lhs_vec(bin_count), rhs_vec(bin_count);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         lhs_vec[i].s[0] = lhs[2 * i];
         lhs_vec[i].s[1] = lhs[2 * i + 1];
         rhs_vec[i].s[0] = rhs[2 * i];
         rhs_vec[i].s[1] = rhs[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem lhs_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                         bytes, lhs_vec.data(), &err);
      cl_mem rhs_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                         bytes, rhs_vec.data(), &err);
      cl_mem out_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, bytes, nullptr, &err);
      if(err != CL_SUCCESS)
        {
         if(lhs_buffer) clReleaseMemObject(lhs_buffer);
         if(rhs_buffer) clReleaseMemObject(rhs_buffer);
         if(out_buffer) clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(convolution_kernel_, 0, sizeof(cl_mem), &lhs_buffer);
      err |= clSetKernelArg(convolution_kernel_, 1, sizeof(cl_mem), &rhs_buffer);
      err |= clSetKernelArg(convolution_kernel_, 2, sizeof(cl_mem), &out_buffer);
      err |= clSetKernelArg(convolution_kernel_, 3, sizeof(double), &norm);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(lhs_buffer);
         clReleaseMemObject(rhs_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = bin_count;
      err = clEnqueueNDRangeKernel(queue_, convolution_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(lhs_buffer);
         clReleaseMemObject(rhs_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<cl_double2> output_vec(bin_count);
      err = clEnqueueReadBuffer(queue_, out_buffer, CL_TRUE, 0, bytes, output_vec.data(), 0, nullptr, nullptr);

      clReleaseMemObject(lhs_buffer);
      clReleaseMemObject(rhs_buffer);
      clReleaseMemObject(out_buffer);

      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      out.resize(bin_count * 2);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         out[2 * i]     = output_vec[i].s[0];
         out[2 * i + 1] = output_vec[i].s[1];
        }
      return JobStatus::Ok;
     }

   JobStatus correlation(const double* lhs,
                         const double* rhs,
                         std::size_t   bin_count,
                         std::vector<double>& out) override
     {
      const size_t bytes = bin_count * sizeof(cl_double2);
      std::vector<cl_double2> lhs_vec(bin_count), rhs_vec(bin_count);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         lhs_vec[i].s[0] = lhs[2 * i];
         lhs_vec[i].s[1] = lhs[2 * i + 1];
         rhs_vec[i].s[0] = rhs[2 * i];
         rhs_vec[i].s[1] = rhs[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem lhs_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                         bytes, lhs_vec.data(), &err);
      cl_mem rhs_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                         bytes, rhs_vec.data(), &err);
      cl_mem out_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, bytes, nullptr, &err);
      if(err != CL_SUCCESS)
        {
         if(lhs_buffer) clReleaseMemObject(lhs_buffer);
         if(rhs_buffer) clReleaseMemObject(rhs_buffer);
         if(out_buffer) clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(correlation_kernel_, 0, sizeof(cl_mem), &lhs_buffer);
      err |= clSetKernelArg(correlation_kernel_, 1, sizeof(cl_mem), &rhs_buffer);
      err |= clSetKernelArg(correlation_kernel_, 2, sizeof(cl_mem), &out_buffer);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(lhs_buffer);
         clReleaseMemObject(rhs_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = bin_count;
      err = clEnqueueNDRangeKernel(queue_, correlation_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(lhs_buffer);
         clReleaseMemObject(rhs_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<cl_double2> output_vec(bin_count);
      err = clEnqueueReadBuffer(queue_, out_buffer, CL_TRUE, 0, bytes, output_vec.data(), 0, nullptr, nullptr);

      clReleaseMemObject(lhs_buffer);
      clReleaseMemObject(rhs_buffer);
      clReleaseMemObject(out_buffer);

      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      out.resize(bin_count * 2);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         out[2 * i]     = output_vec[i].s[0];
         out[2 * i + 1] = output_vec[i].s[1];
        }
      return JobStatus::Ok;
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

      const size_t dst_len = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(static_cast<double>(length) * factor)));
      cl_int err = CL_SUCCESS;
      cl_mem input_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                           length * sizeof(double), const_cast<double*>(input), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, dst_len * sizeof(double), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const int src_len_int = static_cast<int>(length);
      err  = clSetKernelArg(resample_kernel_, 0, sizeof(cl_mem), &input_buffer);
      err |= clSetKernelArg(resample_kernel_, 1, sizeof(cl_mem), &output_buffer);
      err |= clSetKernelArg(resample_kernel_, 2, sizeof(int), &src_len_int);
      err |= clSetKernelArg(resample_kernel_, 3, sizeof(double), &factor);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = dst_len;
      err = clEnqueueNDRangeKernel(queue_, resample_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      out.resize(dst_len);
      err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, dst_len * sizeof(double), out.data(), 0, nullptr, nullptr);

      clReleaseMemObject(input_buffer);
      clReleaseMemObject(output_buffer);

      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      return JobStatus::Ok;
     }

   JobStatus zero_pad(const double* input,
                      std::size_t   length,
                      std::size_t   pad_left,
                      std::size_t   pad_right,
                      std::vector<double>& out) override
     {
      const size_t dst_len = pad_left + length + pad_right;
      cl_int err = CL_SUCCESS;
      cl_mem input_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                           length * sizeof(double), const_cast<double*>(input), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, dst_len * sizeof(double), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const int src_len_int = static_cast<int>(length);
      const int pad_left_int = static_cast<int>(pad_left);
      err  = clSetKernelArg(zero_pad_kernel_, 0, sizeof(cl_mem), &input_buffer);
      err |= clSetKernelArg(zero_pad_kernel_, 1, sizeof(cl_mem), &output_buffer);
      err |= clSetKernelArg(zero_pad_kernel_, 2, sizeof(int), &src_len_int);
      err |= clSetKernelArg(zero_pad_kernel_, 3, sizeof(int), &pad_left_int);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = dst_len;
      err = clEnqueueNDRangeKernel(queue_, zero_pad_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      out.resize(dst_len);
      err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, dst_len * sizeof(double), out.data(), 0, nullptr, nullptr);

      clReleaseMemObject(input_buffer);
      clReleaseMemObject(output_buffer);

      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      return JobStatus::Ok;
     }

   JobStatus remove_dc(const double* input,
                       std::size_t   length,
                       int           mode,
                       double        alpha,
                       std::vector<double>& out) override
     {
      out.resize(length);
      if(length == 0)
         return JobStatus::Ok;

      if(mode == 0)
        {
         double sum = 0.0;
         for(std::size_t i = 0; i < length; ++i)
            sum += input[i];
         const double mean = sum / static_cast<double>(length);

         cl_int err = CL_SUCCESS;
         cl_mem input_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              length * sizeof(double), const_cast<double*>(input), &err);
         if(err != CL_SUCCESS)
            return JobStatus::ErrorBackendFailure;
         cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, length * sizeof(double), nullptr, &err);
         if(err != CL_SUCCESS)
           {
            clReleaseMemObject(input_buffer);
            return JobStatus::ErrorBackendFailure;
           }

         err  = clSetKernelArg(subtract_mean_kernel_, 0, sizeof(cl_mem), &input_buffer);
         err |= clSetKernelArg(subtract_mean_kernel_, 1, sizeof(cl_mem), &output_buffer);
         err |= clSetKernelArg(subtract_mean_kernel_, 2, sizeof(double), &mean);
         if(err != CL_SUCCESS)
           {
            clReleaseMemObject(input_buffer);
            clReleaseMemObject(output_buffer);
            return JobStatus::ErrorBackendFailure;
           }

         const size_t global = length;
         err = clEnqueueNDRangeKernel(queue_, subtract_mean_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
         if(err != CL_SUCCESS)
           {
            clReleaseMemObject(input_buffer);
            clReleaseMemObject(output_buffer);
            return JobStatus::ErrorBackendFailure;
           }

         err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, length * sizeof(double), out.data(), 0, nullptr, nullptr);
         clReleaseMemObject(input_buffer);
         clReleaseMemObject(output_buffer);
      if(err != CL_SUCCESS)
            return JobStatus::ErrorBackendFailure;

         return JobStatus::Ok;
        }

      if(mode == 1)
        {
         const double a = std::clamp(alpha, 0.0, 1.0);
         cl_int err = CL_SUCCESS;
         cl_mem input_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              length * sizeof(double), const_cast<double*>(input), &err);
         if(err != CL_SUCCESS)
            return JobStatus::ErrorBackendFailure;
         cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, length * sizeof(double), nullptr, &err);
         if(err != CL_SUCCESS)
           {
            clReleaseMemObject(input_buffer);
            return JobStatus::ErrorBackendFailure;
           }

         const int len_int = static_cast<int>(length);
         err  = clSetKernelArg(remove_dc_mode1_kernel_, 0, sizeof(cl_mem), &input_buffer);
         err |= clSetKernelArg(remove_dc_mode1_kernel_, 1, sizeof(cl_mem), &output_buffer);
         err |= clSetKernelArg(remove_dc_mode1_kernel_, 2, sizeof(int), &len_int);
         err |= clSetKernelArg(remove_dc_mode1_kernel_, 3, sizeof(double), &a);
         if(err != CL_SUCCESS)
           {
            clReleaseMemObject(input_buffer);
            clReleaseMemObject(output_buffer);
            return JobStatus::ErrorBackendFailure;
           }

         const size_t global = 1;
         err = clEnqueueNDRangeKernel(queue_, remove_dc_mode1_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
         if(err != CL_SUCCESS)
           {
            clReleaseMemObject(input_buffer);
            clReleaseMemObject(output_buffer);
            return JobStatus::ErrorBackendFailure;
           }

         err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, length * sizeof(double), out.data(), 0, nullptr, nullptr);
         clReleaseMemObject(input_buffer);
         clReleaseMemObject(output_buffer);
         if(err != CL_SUCCESS)
            return JobStatus::ErrorBackendFailure;

         return JobStatus::Ok;
        }

      return JobStatus::ErrorInvalidArgs;
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

      std::vector<cl_double2> spectrum_vec(bin_count);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         spectrum_vec[i].s[0] = spectrum[2 * i];
         spectrum_vec[i].s[1] = spectrum[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem spectrum_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              bin_count * sizeof(cl_double2), spectrum_vec.data(), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem phase_buffer = clCreateBuffer(context_, CL_MEM_READ_WRITE, bin_count * sizeof(double), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         return JobStatus::ErrorBackendFailure;
        }
      cl_mem delta_buffer = clCreateBuffer(context_, CL_MEM_READ_WRITE, bin_count * sizeof(cl_int), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         return JobStatus::ErrorBackendFailure;
        }
      cl_mem wrap_buffer = clCreateBuffer(context_, CL_MEM_READ_WRITE, bin_count * sizeof(cl_int), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         return JobStatus::ErrorBackendFailure;
        }
      cl_mem out_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, bin_count * sizeof(double), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const int bins_int = static_cast<int>(bin_count);
      err  = clSetKernelArg(compute_phase_kernel_, 0, sizeof(cl_mem), &spectrum_buffer);
      err |= clSetKernelArg(compute_phase_kernel_, 1, sizeof(cl_mem), &phase_buffer);
      err |= clSetKernelArg(compute_phase_kernel_, 2, sizeof(int), &bins_int);
      if(err == CL_SUCCESS)
        {
         const size_t global = bin_count;
         err = clEnqueueNDRangeKernel(queue_, compute_phase_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
        }
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(wrap_delta_kernel_, 0, sizeof(cl_mem), &phase_buffer);
      err |= clSetKernelArg(wrap_delta_kernel_, 1, sizeof(cl_mem), &delta_buffer);
      err |= clSetKernelArg(wrap_delta_kernel_, 2, sizeof(int), &bins_int);
      if(err == CL_SUCCESS)
        {
         const size_t global = bin_count;
         err = clEnqueueNDRangeKernel(queue_, wrap_delta_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
        }
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t local = std::min<std::size_t>(256, static_cast<std::size_t>(bin_count));
      const size_t global_scan = ((bin_count + local - 1) / local) * local;
      const int    group_count = static_cast<int>((global_scan == 0 || local == 0) ? 0 : (global_scan / local));

      cl_mem block_totals_buffer = clCreateBuffer(context_, CL_MEM_READ_WRITE, std::max(1, group_count) * sizeof(cl_int), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }
      cl_mem block_offsets_buffer = clCreateBuffer(context_, CL_MEM_READ_WRITE, std::max(1, group_count) * sizeof(cl_int), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(block_totals_buffer);
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(exclusive_scan_kernel_, 0, sizeof(cl_mem), &delta_buffer);
      err |= clSetKernelArg(exclusive_scan_kernel_, 1, sizeof(cl_mem), &wrap_buffer);
      err |= clSetKernelArg(exclusive_scan_kernel_, 2, sizeof(cl_mem), &block_totals_buffer);
      err |= clSetKernelArg(exclusive_scan_kernel_, 3, sizeof(int), &bins_int);
      err |= clSetKernelArg(exclusive_scan_kernel_, 4, local * sizeof(cl_int), nullptr);
      if(err == CL_SUCCESS)
        {
         err = clEnqueueNDRangeKernel(queue_, exclusive_scan_kernel_, 1, nullptr, &global_scan, &local, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
        }
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(block_offsets_buffer);
         clReleaseMemObject(block_totals_buffer);
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<int> block_totals(std::max(1, group_count));
      err = clEnqueueReadBuffer(queue_, block_totals_buffer, CL_TRUE, 0,
                                static_cast<size_t>(std::max(1, group_count)) * sizeof(cl_int),
                                block_totals.data(), 0, nullptr, nullptr);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(block_offsets_buffer);
         clReleaseMemObject(block_totals_buffer);
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<int> block_offsets(std::max(1, group_count), 0);
      int running = 0;
      for(int i = 0; i < group_count; ++i)
        {
         block_offsets[i] = running;
         running += block_totals[i];
        }
      err = clEnqueueWriteBuffer(queue_, block_offsets_buffer, CL_TRUE, 0,
                                 static_cast<size_t>(std::max(1, group_count)) * sizeof(cl_int),
                                 block_offsets.data(), 0, nullptr, nullptr);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(block_offsets_buffer);
         clReleaseMemObject(block_totals_buffer);
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(exclusive_scan_add_kernel_, 0, sizeof(cl_mem), &wrap_buffer);
      err |= clSetKernelArg(exclusive_scan_add_kernel_, 1, sizeof(cl_mem), &block_offsets_buffer);
      err |= clSetKernelArg(exclusive_scan_add_kernel_, 2, sizeof(int), &bins_int);
      if(err == CL_SUCCESS)
        {
         err = clEnqueueNDRangeKernel(queue_, exclusive_scan_add_kernel_, 1, nullptr, &global_scan, &local, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
        }
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(block_offsets_buffer);
         clReleaseMemObject(block_totals_buffer);
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(phase_buffer);
         clReleaseMemObject(delta_buffer);
         clReleaseMemObject(wrap_buffer);
         clReleaseMemObject(out_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(finalize_unwrap_kernel_, 0, sizeof(cl_mem), &phase_buffer);
      err |= clSetKernelArg(finalize_unwrap_kernel_, 1, sizeof(cl_mem), &wrap_buffer);
      err |= clSetKernelArg(finalize_unwrap_kernel_, 2, sizeof(cl_mem), &out_buffer);
      err |= clSetKernelArg(finalize_unwrap_kernel_, 3, sizeof(int), &bins_int);
      if(err == CL_SUCCESS)
        {
         const size_t global = bin_count;
         err = clEnqueueNDRangeKernel(queue_, finalize_unwrap_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
        }

      JobStatus status = JobStatus::Ok;
      if(err == CL_SUCCESS)
        {
         out.resize(bin_count);
         err = clEnqueueReadBuffer(queue_, out_buffer, CL_TRUE, 0, bin_count * sizeof(double), out.data(), 0, nullptr, nullptr);
         if(err != CL_SUCCESS)
            status = JobStatus::ErrorBackendFailure;
        }
      else
        status = JobStatus::ErrorBackendFailure;

      clReleaseMemObject(block_offsets_buffer);
      clReleaseMemObject(block_totals_buffer);
      clReleaseMemObject(spectrum_buffer);
      clReleaseMemObject(phase_buffer);
      clReleaseMemObject(delta_buffer);
      clReleaseMemObject(wrap_buffer);
      clReleaseMemObject(out_buffer);

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

      cl_int err = CL_SUCCESS;
      cl_mem series_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                            length * sizeof(double), const_cast<double*>(time_series), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem metrics_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY,
                                             length * sizeof(wave_pipe::InstantMetricsPayload), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(series_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const int len_int = static_cast<int>(length);
      err  = clSetKernelArg(instant_metrics_kernel_, 0, sizeof(cl_mem), &series_buffer);
      err |= clSetKernelArg(instant_metrics_kernel_, 1, sizeof(cl_mem), &metrics_buffer);
      err |= clSetKernelArg(instant_metrics_kernel_, 2, sizeof(int), &len_int);
      err |= clSetKernelArg(instant_metrics_kernel_, 3, sizeof(int), &smooth_window);
      err |= clSetKernelArg(instant_metrics_kernel_, 4, sizeof(double), &epsilon);
      if(err == CL_SUCCESS)
        {
         const size_t global = length;
         err = clEnqueueNDRangeKernel(queue_, instant_metrics_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
        }

      JobStatus status = JobStatus::Ok;
      if(err == CL_SUCCESS)
        {
         instant_out.resize(length);
         err = clEnqueueReadBuffer(queue_, metrics_buffer, CL_TRUE, 0,
                                   length * sizeof(wave_pipe::InstantMetricsPayload),
                                   instant_out.data(), 0, nullptr, nullptr);
         if(err != CL_SUCCESS)
            status = JobStatus::ErrorBackendFailure;
        }
      else
        status = JobStatus::ErrorBackendFailure;

      clReleaseMemObject(series_buffer);
      clReleaseMemObject(metrics_buffer);
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

      cl_int err = CL_SUCCESS;
      cl_mem metrics_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                             length * sizeof(double), const_cast<double*>(metrics), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem transitions_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY,
                                                 length * sizeof(wave_pipe::TransitionPayload), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(metrics_buffer);
         return JobStatus::ErrorBackendFailure;
        }
      cl_mem count_buffer = clCreateBuffer(context_, CL_MEM_READ_WRITE, sizeof(cl_int), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(metrics_buffer);
         clReleaseMemObject(transitions_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const cl_int zero = 0;
      err = clEnqueueWriteBuffer(queue_, count_buffer, CL_TRUE, 0, sizeof(cl_int), &zero, 0, nullptr, nullptr);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(metrics_buffer);
         clReleaseMemObject(transitions_buffer);
         clReleaseMemObject(count_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const int len_int = static_cast<int>(length);
      err  = clSetKernelArg(detect_transitions_kernel_, 0, sizeof(cl_mem), &metrics_buffer);
      err |= clSetKernelArg(detect_transitions_kernel_, 1, sizeof(cl_mem), &transitions_buffer);
      err |= clSetKernelArg(detect_transitions_kernel_, 2, sizeof(cl_mem), &count_buffer);
      err |= clSetKernelArg(detect_transitions_kernel_, 3, sizeof(int), &len_int);
      err |= clSetKernelArg(detect_transitions_kernel_, 4, sizeof(double), &energy_threshold);
      err |= clSetKernelArg(detect_transitions_kernel_, 5, sizeof(double), &phase_jump_threshold);
      err |= clSetKernelArg(detect_transitions_kernel_, 6, sizeof(int), &min_duration);
      err |= clSetKernelArg(detect_transitions_kernel_, 7, sizeof(int), &type_mask);
      if(err == CL_SUCCESS)
        {
         const size_t local = std::min<std::size_t>(128, static_cast<std::size_t>(length > 0 ? length : 1));
         const size_t global = ((length + local - 1) / local) * local;
         err = clEnqueueNDRangeKernel(queue_, detect_transitions_kernel_, 1, nullptr, &global, &local, 0, nullptr, nullptr);
         if(err == CL_SUCCESS) err = clFinish(queue_);
        }

      JobStatus status = JobStatus::Ok;
      cl_int host_count = 0;
      if(err == CL_SUCCESS)
        {
         err = clEnqueueReadBuffer(queue_, count_buffer, CL_TRUE, 0, sizeof(cl_int), &host_count, 0, nullptr, nullptr);
         if(err != CL_SUCCESS)
            status = JobStatus::ErrorBackendFailure;
        }
      else
        status = JobStatus::ErrorBackendFailure;

      host_count = std::max(0, std::min(host_count, static_cast<cl_int>(length)));
      transitions.resize(static_cast<std::size_t>(host_count));
      if(status == JobStatus::Ok && host_count > 0)
        {
         err = clEnqueueReadBuffer(queue_, transitions_buffer, CL_TRUE, 0,
                                   static_cast<std::size_t>(host_count) * sizeof(wave_pipe::TransitionPayload),
                                   transitions.data(), 0, nullptr, nullptr);
         if(err != CL_SUCCESS)
            status = JobStatus::ErrorBackendFailure;
        }

      clReleaseMemObject(metrics_buffer);
      clReleaseMemObject(transitions_buffer);
      clReleaseMemObject(count_buffer);

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

      std::vector<cl_double2> spectrum_vec(bin_count);
      for(std::size_t i = 0; i < bin_count; ++i)
        {
         spectrum_vec[i].s[0] = spectrum[2 * i];
         spectrum_vec[i].s[1] = spectrum[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem spectrum_buffer = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              bin_count * sizeof(cl_double2), spectrum_vec.data(), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;
      cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, bin_count * sizeof(double), nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(hartley_kernel_, 0, sizeof(cl_mem), &spectrum_buffer);
      err |= clSetKernelArg(hartley_kernel_, 1, sizeof(cl_mem), &output_buffer);
      err |= clSetKernelArg(hartley_kernel_, 2, sizeof(double), &scale);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = bin_count;
      err = clEnqueueNDRangeKernel(queue_, hartley_kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS) err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(spectrum_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      out.resize(bin_count);
      err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, bin_count * sizeof(double), out.data(), 0, nullptr, nullptr);
      clReleaseMemObject(spectrum_buffer);
      clReleaseMemObject(output_buffer);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      return JobStatus::Ok;
     }

  private:
   cl_platform_id  platform_{};
   cl_device_id    device_{};
   cl_context      context_{};
   cl_command_queue queue_{};
   cl_program      program_{};

   cl_kernel mask_complex_kernel_{};
   cl_kernel mask_real_kernel_{};
   cl_kernel denoise_kernel_{};
   cl_kernel upscale_kernel_{};
   cl_kernel downscale_kernel_{};
   cl_kernel convolution_kernel_{};
   cl_kernel correlation_kernel_{};
   cl_kernel resample_kernel_{};
   cl_kernel zero_pad_kernel_{};
   cl_kernel subtract_mean_kernel_{};
   cl_kernel hartley_kernel_{};
   cl_kernel remove_dc_mode1_kernel_{};
   cl_kernel compute_phase_kernel_{};
   cl_kernel wrap_delta_kernel_{};
   cl_kernel exclusive_scan_kernel_{};
   cl_kernel exclusive_scan_add_kernel_{};
   cl_kernel finalize_unwrap_kernel_{};
   cl_kernel instant_metrics_kernel_{};
   cl_kernel detect_transitions_kernel_{};
  };
} // namespace

std::unique_ptr<IKernelExecutor> CreateOpenClKernelExecutor(const GpuConfig& cfg)
  {
   try
     {
      return std::make_unique<OpenClKernelExecutor>(cfg);
     }
   catch(...)
     {
      return nullptr;
     }
  }

} // namespace alglib_gpu

#endif // ALGLIB_GPU_ENABLE_OPENCL
