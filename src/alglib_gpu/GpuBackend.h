#pragma once

#include "GpuCommon.h"
#include "alglib_pipe_messages.h"

#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <unordered_map>
#include <vector>

namespace alglib_gpu
{
class IFftExecutor
  {
public:
   virtual ~IFftExecutor() = default;

   virtual JobStatus execute_complex(std::vector<double>& inout, int fft_len, bool inverse) = 0;
   virtual JobStatus execute_real(const std::vector<double>& input,
                                  std::vector<double>& output,
                                  int fft_len,
                                  bool inverse) = 0;
  };

std::unique_ptr<IFftExecutor> CreateCudaExecutor(const GpuConfig& cfg);
std::unique_ptr<IFftExecutor> CreateOpenClExecutor(const GpuConfig& cfg);

class IKernelExecutor
  {
public:
   virtual ~IKernelExecutor() = default;

   virtual JobStatus apply_mask(const double* spectrum,
                                std::size_t   bin_count,
                                const double* mask,
                                std::size_t   mask_len,
                                bool          mask_is_complex,
                                int           mode,
                                std::vector<double>& out) = 0;

   virtual JobStatus denoise(const double* spectrum,
                              std::size_t   bin_count,
                              int           method,
                              double        threshold,
                              double        beta,
                              int           iterations,
                              std::vector<double>& out) = 0;

   virtual JobStatus upscale(const double* spectrum,
                              std::size_t   bin_count,
                              double        factor,
                              int           mode,
                              int           normalize,
                              std::vector<double>& out) = 0;

   virtual JobStatus downscale(const double* spectrum,
                                std::size_t   bin_count,
                                double        factor,
                                int           mode,
                                int           anti_alias,
                                std::vector<double>& out) = 0;

   virtual JobStatus convolution(const double* lhs,
                                  const double* rhs,
                                  std::size_t   bin_count,
                                  int           normalize,
                                  std::vector<double>& out) = 0;

   virtual JobStatus correlation(const double* lhs,
                                  const double* rhs,
                                  std::size_t   bin_count,
                                  std::vector<double>& out) = 0;

   virtual JobStatus resample(const double* input,
                              std::size_t   length,
                              double        factor,
                              double        cutoff,
                              int           method,
                              std::vector<double>& out) = 0;

   virtual JobStatus zero_pad(const double* input,
                              std::size_t   length,
                              std::size_t   pad_left,
                              std::size_t   pad_right,
                              std::vector<double>& out) = 0;

   virtual JobStatus remove_dc(const double* input,
                               std::size_t   length,
                               int           mode,
                               double        alpha,
                               std::vector<double>& out) = 0;

   virtual JobStatus phase_unwrap(const double* spectrum,
                                  std::size_t   bin_count,
                                  int           method,
                                  std::vector<double>& out) = 0;

   virtual JobStatus instant_metrics(const double* time_series,
                                     std::size_t   length,
                                     int           smooth_window,
                                     double        epsilon,
                                     std::vector<wave_pipe::InstantMetricsPayload>& instant_out) = 0;

  virtual JobStatus detect_transitions(const double* metrics,
                                       std::size_t   length,
                                       double        energy_threshold,
                                       double        phase_jump_threshold,
                                       int           min_duration,
                                       int           type_mask,
                                       std::vector<wave_pipe::TransitionPayload>& transitions) = 0;

   virtual JobStatus hartley_from_complex(const double* spectrum,
                                          std::size_t   bin_count,
                                          double        scale,
                                          std::vector<double>& out) = 0;
  };

std::unique_ptr<IKernelExecutor> CreateCudaKernelExecutor(const GpuConfig& cfg);
std::unique_ptr<IKernelExecutor> CreateOpenClKernelExecutor(const GpuConfig& cfg);

class BackendContext
  {
public:
   BackendContext();
   ~BackendContext();

   JobStatus initialize(const GpuConfig& cfg);
   void      shutdown();

   JobStatus execute_fft(const FftSubmitDesc& desc, std::vector<double>& buffer);

   JobStatus apply_mask(const double* spectrum,
                        std::size_t   bin_count,
                        const double* mask,
                        std::size_t   mask_len,
                        bool          mask_is_complex,
                        int           mode,
                        std::vector<double>& out);

   JobStatus spectral_denoise(const double* spectrum,
                              std::size_t   bin_count,
                              int           method,
                              double        threshold,
                              double        beta,
                              int           iterations,
                              std::vector<double>& out);

   JobStatus spectral_upscale(const double* spectrum,
                              std::size_t   bin_count,
                              double        factor,
                              int           mode,
                              int           normalize,
                              std::vector<double>& out);

   JobStatus spectral_downscale(const double* spectrum,
                                std::size_t   bin_count,
                                double        factor,
                                int           mode,
                                int           anti_alias,
                                std::vector<double>& out);

   JobStatus spectral_convolution(const double* lhs,
                                  const double* rhs,
                                  std::size_t   bin_count,
                                  int           normalize,
                                  std::vector<double>& out);

   JobStatus spectral_correlation(const double* lhs,
                                  const double* rhs,
                                  std::size_t   bin_count,
                                  std::vector<double>& out);

   JobStatus resample_time_series(const double* input,
                                  std::size_t   length,
                                  double        factor,
                                  double        cutoff,
                                  int           method,
                                  std::vector<double>& out);

   JobStatus zero_pad_time_series(const double* input,
                                  std::size_t   length,
                                  std::size_t   pad_left,
                                  std::size_t   pad_right,
                                  std::vector<double>& out);

   JobStatus remove_dc_time_series(const double* input,
                                   std::size_t   length,
                                   int           mode,
                                   double        alpha,
                                   std::vector<double>& out);

   JobStatus spectral_phase_unwrap(const double* spectrum,
                                   std::size_t   bin_count,
                                   int           method,
                                   std::vector<double>& out);

   JobStatus spectral_instant_metrics(const double* time_series,
                                      std::size_t   length,
                                      int           smooth_window,
                                      double        epsilon,
                                      std::vector<wave_pipe::InstantMetricsPayload>& instant);

  JobStatus spectral_detect_transitions(const double* metrics,
                                        std::size_t   length,
                                        double        energy_threshold,
                                        double        phase_jump_threshold,
                                        int           min_duration,
                                        int           type_mask,
                                        std::vector<wave_pipe::TransitionPayload>& transitions);

   JobStatus hartley_transform(const double* input,
                               std::size_t   length,
                               bool          inverse,
                               std::vector<double>& out);

  BackendType resolved_backend() const;

 private:
   std::mutex                      mutex_;
   GpuConfig                       config_;
   BackendType                     resolved_backend_;
   std::unique_ptr<IFftExecutor>   executor_;
   std::unique_ptr<IKernelExecutor> kernel_executor_;
   bool                            initialized_;
  };
} // namespace alglib_gpu
