#include "GpuBackend.h"

#ifdef ALGLIB_GPU_ENABLE_CUDA

#include <cuda_runtime.h>
#include <cufft.h>

#include <atomic>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <vector>

inline alglib_gpu::JobStatus translate_cufft(cufftResult res)
  {
   using alglib_gpu::JobStatus;
   return (res == CUFFT_SUCCESS) ? JobStatus::Ok : JobStatus::ErrorBackendFailure;
  }

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


struct PlanKey
  {
   int     fft_len;
   FftKind kind;

   bool operator==(const PlanKey& other) const
     {
      return fft_len == other.fft_len && kind == other.kind;
     }
  };

struct PlanKeyHasher
  {
   std::size_t operator()(const PlanKey& key) const noexcept
     {
      return std::hash<int>{}(key.fft_len) ^ (static_cast<std::size_t>(key.kind) << 1);
     }
  };

class CudaFftExecutor final : public IFftExecutor
  {
public:
   explicit CudaFftExecutor(const GpuConfig& cfg)
       : config_(cfg),
         next_stream_(0)
     {
      const int device = std::max(0, cfg.device_index);
      cudaSetDevice(device);

      const int stream_count = cfg.stream_count > 0 ? cfg.stream_count : 1;
      streams_.reserve(stream_count);
      for(int i = 0; i < stream_count; ++i)
        {
         cudaStream_t stream{};
         if(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) != cudaSuccess)
            throw std::runtime_error("Failed to create CUDA stream");
         streams_.push_back(stream);
        }
     }

   ~CudaFftExecutor() override
     {
      for(auto stream : streams_)
        {
         if(stream != nullptr)
            cudaStreamDestroy(stream);
        }
      for(auto& [key, plan] : plan_cache_)
        cufftDestroy(plan);
     }

   JobStatus execute_complex(std::vector<double>& inout, int fft_len, bool inverse) override
     {
      const size_t bytes = static_cast<size_t>(fft_len) * sizeof(cufftDoubleComplex);

      cufftDoubleComplex* d_data = nullptr;
      auto status                = translate_cuda(cudaMalloc(reinterpret_cast<void**>(&d_data), bytes));
      if(status != JobStatus::Ok)
         return status;

      std::vector<cufftDoubleComplex> host(fft_len);
      for(int i = 0; i < fft_len; ++i)
        {
         host[i].x = inout[2 * i];
         host[i].y = inout[2 * i + 1];
        }

      status = translate_cuda(cudaMemcpy(d_data, host.data(), bytes, cudaMemcpyHostToDevice));
      if(status != JobStatus::Ok)
        {
         cudaFree(d_data);
         return status;
        }

      cufftHandle plan = obtain_plan(fft_len, FftKind::Complex);
      auto        stream = pick_stream();
      cufftSetStream(plan, stream);

      const int direction = inverse ? CUFFT_INVERSE : CUFFT_FORWARD;
      auto      cufft_err = cufftExecZ2Z(plan, d_data, d_data, direction);
      auto      cufft_st  = translate_cufft(cufft_err);
      if(cufft_st != JobStatus::Ok)
        {
         cudaFree(d_data);
         return cufft_st;
        }

      status = translate_cuda(cudaMemcpy(inout.data(), d_data, bytes, cudaMemcpyDeviceToHost));
      cudaFree(d_data);
      if(status != JobStatus::Ok)
         return status;

      // cuFFT does not normalise inverse transform; mimic ALGLIB (divide by fft_len).
      if(inverse)
        {
         const double scale = 1.0 / static_cast<double>(fft_len);
         for(int i = 0; i < fft_len; ++i)
           {
            inout[2 * i] *= scale;
            inout[2 * i + 1] *= scale;
           }
        }

      return JobStatus::Ok;
     }

   JobStatus execute_real(const std::vector<double>& input,
                          std::vector<double>&       output,
                          int                        fft_len,
                          bool                       inverse) override
     {
      if(fft_len <= 0)
         return JobStatus::ErrorInvalidArgs;

      if(!inverse)
        {
         if(static_cast<int>(input.size()) != fft_len)
            return JobStatus::ErrorInvalidArgs;

         std::vector<double> complex_data(2 * fft_len, 0.0);
         for(int i = 0; i < fft_len; ++i)
            complex_data[2 * i] = input[i];

         auto status = execute_complex(complex_data, fft_len, false);
         if(status != JobStatus::Ok)
            return status;

         output = std::move(complex_data);
         return JobStatus::Ok;
        }

      if(static_cast<int>(input.size()) != 2 * fft_len)
         return JobStatus::ErrorInvalidArgs;

      std::vector<double> complex_data = input;
      auto                 status       = execute_complex(complex_data, fft_len, true);
      if(status != JobStatus::Ok)
         return status;

      output.resize(fft_len);
      for(int i = 0; i < fft_len; ++i)
         output[i] = complex_data[2 * i];

      return JobStatus::Ok;
     }

private:
   cufftHandle obtain_plan(int fft_len, FftKind kind)
     {
      PlanKey key{fft_len, kind};
      std::lock_guard<std::mutex> lock(plan_mutex_);
      auto                        it = plan_cache_.find(key);
      if(it != plan_cache_.end())
         return it->second;

      cufftHandle plan{};
      cufftResult res = (kind == FftKind::Complex)
                            ? cufftPlan1d(&plan, fft_len, CUFFT_Z2Z, 1)
                            : cufftPlan1d(&plan, fft_len, CUFFT_D2Z, 1);
      if(res != CUFFT_SUCCESS)
         throw std::runtime_error("Failed to create cuFFT plan");

      plan_cache_.emplace(key, plan);
      return plan;
     }

   cudaStream_t pick_stream()
     {
      if(streams_.empty())
         return nullptr;

      const auto index = next_stream_.fetch_add(1, std::memory_order_relaxed) % streams_.size();
      return streams_[index];
     }

private:
   GpuConfig                                         config_;
   std::vector<cudaStream_t>                         streams_;
   std::unordered_map<PlanKey, cufftHandle, PlanKeyHasher> plan_cache_;
   std::mutex                                        plan_mutex_;
   std::atomic<size_t>                               next_stream_;
  };
} // namespace

std::unique_ptr<IFftExecutor> CreateCudaExecutor(const GpuConfig& cfg)
  {
   try
     {
      return std::make_unique<CudaFftExecutor>(cfg);
     }
   catch(...)
     {
      return nullptr;
     }
  }
} // namespace alglib_gpu

#endif // ALGLIB_GPU_ENABLE_CUDA
