#pragma once

#include <cstddef>
#include <cstdint>

#ifdef _WIN32
  #ifdef ALGLIB_GPU_BUILD
    #define ALGLIB_GPU_API extern "C" __declspec(dllexport)
  #else
    #define ALGLIB_GPU_API extern "C" __declspec(dllimport)
  #endif
#else
  #define ALGLIB_GPU_API extern "C"
#endif

namespace alglib_gpu
{
enum class BackendType : int
  {
   Auto = 0,
   CpuFallback = 1,
   Cuda = 2,
   OpenCL = 3
  };

enum class FftKind : int
  {
   Complex = 0,
   Real = 1
  };

enum class JobStatus : int
  {
   Ok = 0,
   Pending = 1,
   ErrorGeneric = -1,
   ErrorInvalidArgs = -2,
   ErrorNoBackend = -3,
   ErrorBackendFailure = -4,
   ErrorTimeout = -5
  };

struct alignas(8) GpuConfig
  {
   BackendType backend;
   int         device_index;
   int         stream_count;
   int         min_gpu_window;
   int         job_timeout_ms;
  };

struct FftSubmitDesc
  {
   const double* data;
   int           fft_len;
   FftKind       kind;
   bool          inverse;
  };
} // namespace alglib_gpu

struct AlglibGpuConfig
  {
   int backend;          // cast to BackendType
   int device_index;
   int stream_count;
   int min_gpu_window;
   int job_timeout_ms;
};

struct AlglibGpuJobHandle
  {
   std::uint64_t value;
};
