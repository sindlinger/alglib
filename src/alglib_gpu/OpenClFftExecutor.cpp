#include "GpuBackend.h"

#ifdef ALGLIB_GPU_ENABLE_OPENCL

// Suprimir deprecations do OpenCL 1.2 para compilar limpo com /WX
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
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace alglib_gpu
{
namespace
{
constexpr const char* kKernelSource = R"(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void dft(__global const double2* input,
                  __global double2*       output,
                  const int               n,
                  const int               inverse)
{
   const int k = get_global_id(0);
   if(k >= n)
      return;

   const double sign = inverse ? 1.0 : -1.0;
   double2 sum = (double2)(0.0, 0.0);

   for(int t = 0; t < n; ++t)
     {
      const double angle = sign * 6.283185307179586476925286766559 * (double)(k * t) / (double)n;
      const double c = cos(angle);
      const double s = sin(angle);
      const double2 val = input[t];
      sum.x += val.x * c - val.y * s;
      sum.y += val.x * s + val.y * c;
     }

   if(inverse)
     {
      sum.x /= (double)n;
      sum.y /= (double)n;
     }

   output[k] = sum;
}
)";

class OpenClError : public std::runtime_error
  {
public:
   explicit OpenClError(const std::string& msg, cl_int code)
       : std::runtime_error(msg + " (err = " + std::to_string(code) + ")")
     {
     }
  };

class OpenClFftExecutor final : public IFftExecutor
  {
public:
   explicit OpenClFftExecutor(const GpuConfig& cfg)
     {
      cl_uint platform_count = 0;
      cl_int  err            = clGetPlatformIDs(0, nullptr, &platform_count);
      if(err != CL_SUCCESS || platform_count == 0)
         throw OpenClError("Nenhuma plataforma OpenCL disponível", err);

      std::vector<cl_platform_id> platforms(platform_count);
      err = clGetPlatformIDs(platform_count, platforms.data(), nullptr);
      if(err != CL_SUCCESS)
         throw OpenClError("Falha ao listar plataformas OpenCL", err);

      struct DeviceInfo
        {
         cl_platform_id platform;
         cl_device_id   device;
        };

      std::vector<DeviceInfo> devices;
      for(auto platform : platforms)
        {
         cl_uint device_count = 0;
         err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &device_count);
         if(err != CL_SUCCESS || device_count == 0)
            continue;

         std::vector<cl_device_id> local(device_count);
         err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, device_count, local.data(), nullptr);
         if(err == CL_SUCCESS)
           {
            for(auto device : local)
               devices.push_back(DeviceInfo{platform, device});
           }
        }

      if(devices.empty())
        throw std::runtime_error("Nenhum dispositivo OpenCL GPU encontrado para fallback.");

      const size_t index = std::min(static_cast<size_t>(std::max(0, cfg.device_index)),
                                    devices.size() - 1);
      platform_ = devices[index].platform;
      device_   = devices[index].device;

      context_ = clCreateContext(nullptr, 1, &device_, nullptr, nullptr, &err);
      if(err != CL_SUCCESS)
         throw OpenClError("Falha ao criar contexto OpenCL", err);

      queue_ = clCreateCommandQueue(context_, device_, 0, &err);
      if(err != CL_SUCCESS)
         throw OpenClError("Falha ao criar command queue OpenCL", err);

      const char* source_ptr = kKernelSource;
      const size_t source_len = std::strlen(kKernelSource);
      program_ = clCreateProgramWithSource(context_, 1, &source_ptr, &source_len, &err);
      if(err != CL_SUCCESS)
         throw OpenClError("Falha ao criar programa OpenCL", err);

      err = clBuildProgram(program_, 1, &device_, nullptr, nullptr, nullptr);
      if(err != CL_SUCCESS)
        {
         size_t log_size = 0;
         clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
         std::string log(log_size, '\0');
         clGetProgramBuildInfo(program_, device_, CL_PROGRAM_BUILD_LOG, log_size, log.data(), nullptr);
         throw std::runtime_error("Erro ao compilar kernel OpenCL:\n" + log);
        }

      kernel_ = clCreateKernel(program_, "dft", &err);
      if(err != CL_SUCCESS)
         throw OpenClError("Falha ao criar kernel DFT OpenCL", err);
     }

   ~OpenClFftExecutor() override
     {
      if(kernel_)
         clReleaseKernel(kernel_);
      if(program_)
         clReleaseProgram(program_);
      if(queue_)
         clReleaseCommandQueue(queue_);
      if(context_)
         clReleaseContext(context_);
     }

   JobStatus execute_complex(std::vector<double>& inout, int fft_len, bool inverse) override
     {
      if(static_cast<int>(inout.size()) != 2 * fft_len)
         return JobStatus::ErrorInvalidArgs;

      const size_t bytes = static_cast<size_t>(fft_len) * sizeof(cl_double2);
      std::vector<cl_double2> host_in(fft_len);
      for(int i = 0; i < fft_len; ++i)
        {
         host_in[i].s[0] = inout[2 * i];
         host_in[i].s[1] = inout[2 * i + 1];
        }

      cl_int err = CL_SUCCESS;
      cl_mem input_buffer =
          clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, bytes, host_in.data(), &err);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      cl_mem output_buffer = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, bytes, nullptr, &err);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      err  = clSetKernelArg(kernel_, 0, sizeof(cl_mem), &input_buffer);
      err |= clSetKernelArg(kernel_, 1, sizeof(cl_mem), &output_buffer);
      err |= clSetKernelArg(kernel_, 2, sizeof(int), &fft_len);
      const int inverse_flag = inverse ? 1 : 0;
      err |= clSetKernelArg(kernel_, 3, sizeof(int), &inverse_flag);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      const size_t global = static_cast<size_t>(fft_len);
      err = clEnqueueNDRangeKernel(queue_, kernel_, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
      if(err == CL_SUCCESS)
         err = clFinish(queue_);
      if(err != CL_SUCCESS)
        {
         clReleaseMemObject(input_buffer);
         clReleaseMemObject(output_buffer);
         return JobStatus::ErrorBackendFailure;
        }

      std::vector<cl_double2> host_out(fft_len);
      err = clEnqueueReadBuffer(queue_, output_buffer, CL_TRUE, 0, bytes, host_out.data(), 0, nullptr, nullptr);
      clReleaseMemObject(input_buffer);
      clReleaseMemObject(output_buffer);
      if(err != CL_SUCCESS)
         return JobStatus::ErrorBackendFailure;

      for(int i = 0; i < fft_len; ++i)
        {
         inout[2 * i]     = host_out[i].s[0];
         inout[2 * i + 1] = host_out[i].s[1];
        }

      return JobStatus::Ok;
     }

   JobStatus execute_real(const std::vector<double>& input,
                          std::vector<double>&       output,
                          int                        fft_len,
                          bool                       inverse) override
     {
      if(!inverse)
        {
         if(static_cast<int>(input.size()) != fft_len)
            return JobStatus::ErrorInvalidArgs;

         std::vector<double> complex_buffer(2 * fft_len, 0.0);
         for(int i = 0; i < fft_len; ++i)
            complex_buffer[2 * i] = input[i];

         auto status = execute_complex(complex_buffer, fft_len, false);
         if(status != JobStatus::Ok)
            return status;

         output = std::move(complex_buffer);
         return JobStatus::Ok;
        }

      if(static_cast<int>(input.size()) != 2 * fft_len)
         return JobStatus::ErrorInvalidArgs;

      std::vector<double> complex_buffer = input;
      auto                 status         = execute_complex(complex_buffer, fft_len, true);
      if(status != JobStatus::Ok)
         return status;

      output.resize(fft_len);
      for(int i = 0; i < fft_len; ++i)
         output[i] = complex_buffer[2 * i];

      return JobStatus::Ok;
     }

private:
   cl_platform_id platform_ = nullptr;
   cl_device_id   device_   = nullptr;
   cl_context     context_  = nullptr;
   cl_command_queue queue_  = nullptr;
   cl_program     program_  = nullptr;
   cl_kernel      kernel_   = nullptr;
  };
} // namespace

std::unique_ptr<IFftExecutor> CreateOpenClExecutor(const GpuConfig& cfg)
  {
   try
     {
      return std::make_unique<OpenClFftExecutor>(cfg);
     }
   catch(...)
     {
      return nullptr;
     }
  }
} // namespace alglib_gpu

#else

namespace alglib_gpu
{
std::unique_ptr<IFftExecutor> CreateOpenClExecutor(const GpuConfig&)
  {
   return nullptr;
  }
} // namespace alglib_gpu

#endif
