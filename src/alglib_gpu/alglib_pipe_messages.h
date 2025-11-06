#pragma once

#include <cstdint>

namespace wave_pipe
{
constexpr std::uint32_t MESSAGE_MAGIC   = 0x4C574650; // 'LWFP'
constexpr std::uint32_t PROTOCOL_VERSION = 1;

enum class Command : std::uint32_t
  {
   PROCESS_SUBMIT   = 20,
  PROCESS_FETCH    = 21,

   BATCH_SUBMIT     = 30,
   BATCH_FETCH      = 31
  };

enum class Operation : std::uint32_t
  {
   PING                  = 0,
   FFT_COMPLEX_FORWARD   = 1,
   FFT_COMPLEX_INVERSE   = 2,
   FFT_REAL_FORWARD      = 3,
   FFT_REAL_INVERSE      = 4,
   FFT_SINE_TRANSFORM    = 5,
   FFT_COSINE_TRANSFORM  = 6,
   FFT_TWO_REAL          = 7,

   SPECTRAL_RESAMPLE     = 10,
   SPECTRAL_ZERO_PAD     = 11,
   SPECTRAL_REMOVE_DC    = 12,
   SPECTRAL_FILTER_GAUSSIAN = 13,
   SPECTRAL_FILTER_NOTCH = 14,
   SPECTRAL_APPLY_MASK   = 15,
   SPECTRAL_DENOISE      = 16,
   SPECTRAL_UPSCALE      = 17,
   SPECTRAL_DOWNSCALE    = 18,
   SPECTRAL_CONVOLUTION  = 19,
   SPECTRAL_CORRELATION  = 20,
   SPECTRAL_ANALYZE      = 21,
   SPECTRAL_DETECT_TRANSITIONS = 22,
   SPECTRAL_INSTANT_METRICS    = 23,
   SPECTRAL_PHASE_UNWRAP  = 24,

   SSA_CREATE             = 200,
   SSA_DESTROY            = 201,
   SSA_CONFIGURE          = 202,
   SSA_CLEAR              = 203,
   SSA_ADD_SEQUENCE       = 204,
   SSA_APPEND_POINT       = 205,
   SSA_ANALYZE_LAST       = 206,
   SSA_FORECAST_LAST      = 207
  };

enum class Status : std::int32_t
  {
   OK      = 0,
   PENDING = 1,
   ERROR   = -1,
   INVALID = -2,
   TIMEOUT = -3
  };

#pragma pack(push, 1)
struct FetchRequest
  {
   std::uint32_t magic;
   std::uint32_t version;
   std::uint32_t cmd;
   std::uint64_t user_tag;
  };

// FftRequest/FftResponse removidos – usar apenas ProcessRequest/ProcessResponse

struct FftCyclePayload
  {
   std::int32_t index;
   double       amplitude;
   double       phase;
   double       frequency;
   double       period;
  };

struct InstantMetricsPayload
  {
   double amplitude;
   double frequency;
   double energy;
   double phase;
   double envelope_up;
   double envelope_down;
  };

struct TransitionPayload
  {
   std::int32_t start_index;
   std::int32_t end_index;
   double       energy_ratio;
   double       phase_shift;
   std::int32_t type;
   std::int32_t reserved;
  };

struct ProcessRequest
  {
   std::uint32_t magic;
   std::uint32_t version;
   std::uint32_t cmd;
   std::uint64_t user_tag;
   std::uint32_t operation;
   std::uint32_t flags;
   std::int32_t  window_len;
   std::int32_t  aux_len;
   std::uint32_t param_size;
  };

struct ProcessResponse
  {
   std::uint32_t magic;
   std::uint32_t version;
   std::uint32_t cmd;
   std::uint64_t user_tag;
   std::int32_t  status;
   std::uint32_t primary_count;
   std::uint32_t secondary_count;
   std::uint32_t cycle_count;
   std::uint32_t instant_count;
   std::uint32_t transition_count;
   std::uint32_t payload_size;
  };

struct BatchHeader
  {
   std::uint32_t magic;
   std::uint32_t version;
   std::uint32_t cmd;
   std::uint64_t user_tag;
   std::uint32_t task_count;
   std::uint32_t flags;
   std::uint64_t data_size;
   std::uint64_t param_size;
  };

struct BatchTaskDescriptor
  {
   std::uint32_t operation;
   std::uint32_t flags;
   std::int32_t  window_len;
   std::int32_t  aux_len;
   std::uint64_t data_offset;
   std::uint64_t param_offset;
   std::uint32_t payload_hint;
   std::uint32_t reserved;
  };

struct BatchResponse
  {
   std::uint32_t magic;
   std::uint32_t version;
   std::uint32_t cmd;
   std::uint64_t user_tag;
   std::int32_t  status;
   std::uint32_t result_count;
   std::uint64_t payload_size;
   std::uint32_t reserved;
  };
#pragma pack(pop)
} // namespace wave_pipe
