// alglib_bridge.cpp — DLL ponte: comunica com o serviço via Named Pipe
#include <windows.h>
#include <cstdint>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>

#ifdef ERROR
#undef ERROR
#endif
#include "alglib_pipe_messages.h"

namespace bridge {
static const wchar_t* kPipeName = L"\\\\.\\pipe\\alglib-wave_pipe";

static bool WriteExact(HANDLE h, const void* data, DWORD bytes) {
  const std::uint8_t* p = static_cast<const std::uint8_t*>(data);
  DWORD off = 0;
  while (off < bytes) {
    DWORD w = 0;
    if (!::WriteFile(h, p + off, bytes - off, &w, nullptr) || w == 0) return false;
    off += w;
  }
  return true;
}

static bool ReadExact(HANDLE h, void* data, DWORD bytes) {
  std::uint8_t* p = static_cast<std::uint8_t*>(data);
  DWORD off = 0;
  while (off < bytes) {
    DWORD r = 0;
    if (!::ReadFile(h, p + off, bytes - off, &r, nullptr) || r == 0) return false;
    off += r;
  }
  return true;
}

static HANDLE OpenPipe(int timeout_ms) {
  // Espera o pipe ficar disponível dentro do timeout
  auto start = std::chrono::steady_clock::now();
  while (true) {
    if (::WaitNamedPipeW(kPipeName, 200)) {
      HANDLE h = ::CreateFileW(kPipeName,
                               GENERIC_READ | GENERIC_WRITE,
                               FILE_SHARE_READ | FILE_SHARE_WRITE,
                               nullptr,
                               OPEN_EXISTING,
                               0,
                               nullptr);
      if (h != INVALID_HANDLE_VALUE) return h;
    }
    if (timeout_ms >= 0) {
      auto now = std::chrono::steady_clock::now();
      auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
      if (ms > timeout_ms) return INVALID_HANDLE_VALUE;
    }
    ::Sleep(25);
  }
}

static std::uint64_t NextTag() {
  static std::uint64_t tag = 1;
  return tag++;
}

static int SubmitFetch(HANDLE h,
                       wave_pipe::Operation op,
                       const double* primary, int primary_len,
                       const double* secondary, int secondary_len,
                       const std::uint8_t* params, std::uint32_t param_len,
                       std::vector<std::uint8_t>& out_payload,
                       std::uint32_t& out_primary_count,
                       std::uint32_t& out_secondary_count,
                       int timeout_ms) {
  using namespace wave_pipe;
  FetchRequest fetch{};
  ProcessRequest req{};
  const std::uint64_t tag = NextTag();
  req.magic = MESSAGE_MAGIC;
  req.version = PROTOCOL_VERSION;
  req.cmd = static_cast<std::uint32_t>(Command::PROCESS_SUBMIT);
  req.user_tag = tag;
  req.operation = static_cast<std::uint32_t>(op);
  req.flags = 0;
  req.window_len = primary_len;
  req.aux_len = secondary_len;
  req.param_size = param_len;

  if (!WriteExact(h, &req, sizeof(req))) return -5; // IO_ERROR
  if (param_len > 0 && params) {
    if (!WriteExact(h, params, param_len)) return -5;
  }
  if (primary_len > 0 && primary) {
    if (!WriteExact(h, primary, static_cast<DWORD>(primary_len * sizeof(double)))) return -5;
  }
  if (secondary_len > 0 && secondary) {
    if (!WriteExact(h, secondary, static_cast<DWORD>(secondary_len * sizeof(double)))) return -5;
  }
  ::FlushFileBuffers(h);

  // Loop até OK/erro, tratando PENDING com FETCH
  auto start = std::chrono::steady_clock::now();
  while (true) {
    ProcessResponse resp{};
    if (!ReadExact(h, &resp, sizeof(resp))) return -5;
    if (resp.magic != MESSAGE_MAGIC || resp.version != PROTOCOL_VERSION || resp.user_tag != tag) return -1; // ERROR

    if (resp.status == static_cast<std::int32_t>(Status::PENDING)) {
      // Timeout?
      if (timeout_ms >= 0) {
        auto now = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        if (ms > timeout_ms) return -3; // TIMEOUT
      }
      FetchRequest fr{};
      fr.magic = MESSAGE_MAGIC; fr.version = PROTOCOL_VERSION; fr.cmd = static_cast<std::uint32_t>(Command::PROCESS_FETCH); fr.user_tag = tag;
      if (!WriteExact(h, &fr, sizeof(fr))) return -5;
      ::FlushFileBuffers(h);
      ::Sleep(1);
      continue;
    }

    if (resp.status != static_cast<std::int32_t>(Status::OK)) return resp.status; // já é código negativo/erro

    out_primary_count = resp.primary_count;
    out_secondary_count = resp.secondary_count;
    out_payload.clear();
    if (resp.payload_size > 0) {
      out_payload.resize(resp.payload_size);
      if (!ReadExact(h, out_payload.data(), resp.payload_size)) return -5;
    }
    return 0;
  }
}
} // namespace bridge

extern "C" {

__declspec(dllexport)
int AlglibBridge_Ping(unsigned char* out_buf, int out_capacity, int* out_len, int timeout_ms) {
  using namespace bridge;
  HANDLE h = OpenPipe(timeout_ms);
  if (h == INVALID_HANDLE_VALUE) return -5; // IO/timeout de abertura
  std::vector<std::uint8_t> payload;
  std::uint32_t pc=0, sc=0;
  int st = SubmitFetch(h, wave_pipe::Operation::PING, nullptr, 0, nullptr, 0, nullptr, 0, payload, pc, sc, timeout_ms);
  ::CloseHandle(h);
  if (st != 0) return st;
  if (out_len) *out_len = static_cast<int>(payload.size());
  if (out_buf && out_capacity > 0 && !payload.empty()) {
    const int n = static_cast<int>(payload.size());
    const int c = (n < out_capacity) ? n : out_capacity;
    std::memcpy(out_buf, payload.data(), c);
  }
  return 0;
}

// Genérico: envia operação e retorna primary/secondary como doubles;
// o restante do payload (cycles/instant/transitions/metadata) retorna bruto em out_extra.
__declspec(dllexport)
int AlglibBridge_Process(unsigned int operation,
                         const double* primary, int primary_len,
                         const double* secondary, int secondary_len,
                         const unsigned char* params, int param_len,
                         double* out_primary, int out_primary_cap, int* out_primary_len,
                         double* out_secondary, int out_secondary_cap, int* out_secondary_len,
                         unsigned char* out_extra, int out_extra_cap, int* out_extra_len,
                         int timeout_ms) {
  using namespace bridge;
  HANDLE h = OpenPipe(timeout_ms);
  if (h == INVALID_HANDLE_VALUE) return -5;
  std::vector<std::uint8_t> payload;
  std::uint32_t pc=0, sc=0;
  int st = SubmitFetch(h, static_cast<wave_pipe::Operation>(operation),
                       primary, primary_len,
                       secondary, secondary_len,
                       params, static_cast<std::uint32_t>(param_len),
                       payload, pc, sc, timeout_ms);
  ::CloseHandle(h);
  if (st != 0) return st;

  // payload = primary(double[pc]) + secondary(double[sc]) + resto
  const std::size_t bytes_primary   = static_cast<std::size_t>(pc) * sizeof(double);
  const std::size_t bytes_secondary = static_cast<std::size_t>(sc) * sizeof(double);
  std::size_t off = 0;

  if (out_primary_len) *out_primary_len = static_cast<int>(pc);
  if (out_secondary_len) *out_secondary_len = static_cast<int>(sc);
  if (out_extra_len) *out_extra_len = 0;

  if (pc > 0 && out_primary && out_primary_cap > 0) {
    const int to_copy = (pc < static_cast<std::uint32_t>(out_primary_cap)) ? static_cast<int>(pc) : out_primary_cap;
    if (bytes_primary <= payload.size()) std::memcpy(out_primary, payload.data(), static_cast<std::size_t>(to_copy) * sizeof(double));
  }
  off += bytes_primary;
  if (sc > 0 && out_secondary && out_secondary_cap > 0) {
    const int to_copy = (sc < static_cast<std::uint32_t>(out_secondary_cap)) ? static_cast<int>(sc) : out_secondary_cap;
    if (off + bytes_secondary <= payload.size()) std::memcpy(out_secondary, payload.data() + off, static_cast<std::size_t>(to_copy) * sizeof(double));
  }
  off += bytes_secondary;

  if (off < payload.size() && out_extra && out_extra_cap > 0) {
    const int rest = static_cast<int>(payload.size() - off);
    const int to_copy = (rest < out_extra_cap) ? rest : out_extra_cap;
    std::memcpy(out_extra, payload.data() + off, to_copy);
    if (out_extra_len) *out_extra_len = rest;
  } else {
    if (out_extra_len) *out_extra_len = static_cast<int>(payload.size() - off);
  }
  return 0;
}

__declspec(dllexport)
int AlglibBridge_ProcessSubmit(unsigned int operation,
                               const double* primary, int primary_len,
                               const double* secondary, int secondary_len,
                               const unsigned char* params, int param_len,
                               unsigned long long* out_handle,
                               int timeout_ms) {
  using namespace bridge;
  using namespace wave_pipe;
  HANDLE h = OpenPipe(timeout_ms);
  if (h == INVALID_HANDLE_VALUE) return -5;

  const std::uint64_t tag = NextTag();
  ProcessRequest req{};
  req.magic = MESSAGE_MAGIC; req.version = PROTOCOL_VERSION; req.cmd = static_cast<std::uint32_t>(Command::PROCESS_SUBMIT);
  req.user_tag = tag; req.operation = operation; req.flags = 0; req.window_len = primary_len; req.aux_len = secondary_len; req.param_size = static_cast<std::uint32_t>(param_len);

  if (!WriteExact(h, &req, sizeof(req))) { ::CloseHandle(h); return -5; }
  if (param_len > 0 && params) { if (!WriteExact(h, params, static_cast<DWORD>(param_len))) { ::CloseHandle(h); return -5; } }
  if (primary_len > 0 && primary) { if (!WriteExact(h, primary, static_cast<DWORD>(primary_len * sizeof(double)))) { ::CloseHandle(h); return -5; } }
  if (secondary_len > 0 && secondary) { if (!WriteExact(h, secondary, static_cast<DWORD>(secondary_len * sizeof(double)))) { ::CloseHandle(h); return -5; } }
  ::FlushFileBuffers(h);

  // Ler resposta imediata (espera-se PENDING)
  ProcessResponse resp{};
  if (!ReadExact(h, &resp, sizeof(resp))) { ::CloseHandle(h); return -5; }
  ::CloseHandle(h);
  if (resp.magic != MESSAGE_MAGIC || resp.version != PROTOCOL_VERSION || resp.user_tag != tag) return -1;
  if (out_handle) *out_handle = tag;
  // status PENDING esperado; se vier erro, repassar
  if (resp.status == static_cast<std::int32_t>(Status::PENDING)) return 0;
  return resp.status;
}

__declspec(dllexport)
int AlglibBridge_ProcessFetch(unsigned long long handle,
                              double* out_primary, int out_primary_cap, int* out_primary_len,
                              double* out_secondary, int out_secondary_cap, int* out_secondary_len,
                              unsigned char* out_extra, int out_extra_cap, int* out_extra_len,
                              int timeout_ms) {
  using namespace bridge;
  using namespace wave_pipe;
  HANDLE h = OpenPipe(timeout_ms);
  if (h == INVALID_HANDLE_VALUE) return -5;
  auto start = std::chrono::steady_clock::now();
  while (true) {
    FetchRequest fr{}; fr.magic = MESSAGE_MAGIC; fr.version = PROTOCOL_VERSION; fr.cmd = static_cast<std::uint32_t>(Command::PROCESS_FETCH); fr.user_tag = handle;
    if (!WriteExact(h, &fr, sizeof(fr))) { ::CloseHandle(h); return -5; }
    ::FlushFileBuffers(h);

    ProcessResponse resp{};
    if (!ReadExact(h, &resp, sizeof(resp))) { ::CloseHandle(h); return -5; }
    if (resp.magic != MESSAGE_MAGIC || resp.version != PROTOCOL_VERSION || resp.user_tag != handle) { ::CloseHandle(h); return -1; }
    if (resp.status == static_cast<std::int32_t>(Status::PENDING)) {
      if (timeout_ms >= 0) {
        auto now = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        if (ms > timeout_ms) { ::CloseHandle(h); return -3; }
      }
      ::Sleep(2);
      continue;
    }

    // OK ou erro; se erro, retorna status
    if (resp.status != static_cast<std::int32_t>(Status::OK)) { ::CloseHandle(h); return resp.status; }

    // Ler payload
    const std::size_t bytes_primary   = static_cast<std::size_t>(resp.primary_count)   * sizeof(double);
    const std::size_t bytes_secondary = static_cast<std::size_t>(resp.secondary_count) * sizeof(double);
    const std::size_t psize = resp.payload_size;
    std::vector<unsigned char> payload(psize);
    if (psize > 0) {
      if (!ReadExact(h, payload.data(), static_cast<DWORD>(psize))) { ::CloseHandle(h); return -5; }
    }
    ::CloseHandle(h);

    if (out_primary_len) *out_primary_len = static_cast<int>(resp.primary_count);
    if (out_secondary_len) *out_secondary_len = static_cast<int>(resp.secondary_count);
    std::size_t off = 0;
    if (resp.primary_count > 0 && out_primary && out_primary_cap > 0) {
      const int to_copy = (resp.primary_count < static_cast<std::uint32_t>(out_primary_cap)) ? static_cast<int>(resp.primary_count) : out_primary_cap;
      if (bytes_primary <= payload.size()) std::memcpy(out_primary, payload.data(), static_cast<std::size_t>(to_copy) * sizeof(double));
    }
    off += bytes_primary;
    if (resp.secondary_count > 0 && out_secondary && out_secondary_cap > 0) {
      const int to_copy = (resp.secondary_count < static_cast<std::uint32_t>(out_secondary_cap)) ? static_cast<int>(resp.secondary_count) : out_secondary_cap;
      if (off + bytes_secondary <= payload.size()) std::memcpy(out_secondary, payload.data() + off, static_cast<std::size_t>(to_copy) * sizeof(double));
    }
    off += bytes_secondary;
    if (out_extra_len) *out_extra_len = static_cast<int>(psize - off);
    if (off < psize && out_extra && out_extra_cap > 0) {
      const int rest = static_cast<int>(psize - off);
      const int to_copy = (rest < out_extra_cap) ? rest : out_extra_cap;
      std::memcpy(out_extra, payload.data() + off, to_copy);
    }
    return 0;
  }
}

// -------- Batch API --------
__declspec(dllexport)
int AlglibBridge_BatchSubmit(const unsigned int* ops,
                             const unsigned int* task_flags,
                             const int* win_len,
                             const int* aux_len,
                             const unsigned long long* data_off,
                             const unsigned long long* param_off,
                             const unsigned int* payload_hint,
                             int task_count,
                             const double* data_blob, int data_blob_len,
                             const unsigned char* params_blob, int params_blob_len,
                             unsigned long long* out_handle,
                             int timeout_ms) {
  using namespace bridge;
  using namespace wave_pipe;
  HANDLE h = OpenPipe(timeout_ms);
  if (h == INVALID_HANDLE_VALUE) return -5;

  const std::uint64_t tag = NextTag();
  // Send FetchRequest first part of BatchHeader
  FetchRequest fr{}; fr.magic = MESSAGE_MAGIC; fr.version = PROTOCOL_VERSION; fr.cmd = static_cast<std::uint32_t>(Command::BATCH_SUBMIT); fr.user_tag = tag;
  if (!WriteExact(h, &fr, sizeof(fr))) { ::CloseHandle(h); return -5; }

  // Write rest of BatchHeader (task_count, flags, sizes)
  std::uint32_t tcount = static_cast<std::uint32_t>(task_count);
  std::uint32_t flags = 0;
  std::uint64_t dsize = static_cast<std::uint64_t>(data_blob_len * sizeof(double));
  std::uint64_t psize = static_cast<std::uint64_t>(params_blob_len);
  if (!WriteExact(h, &tcount, sizeof(tcount))) { ::CloseHandle(h); return -5; }
  if (!WriteExact(h, &flags, sizeof(flags)))   { ::CloseHandle(h); return -5; }
  if (!WriteExact(h, &dsize, sizeof(dsize)))   { ::CloseHandle(h); return -5; }
  if (!WriteExact(h, &psize, sizeof(psize)))   { ::CloseHandle(h); return -5; }

  // Descriptors
  if (task_count > 0) {
    std::vector<wave_pipe::BatchTaskDescriptor> descs(task_count);
    for (int i=0;i<task_count;++i) {
      descs[i].operation    = ops ? ops[i] : 0u;
      descs[i].flags        = task_flags ? task_flags[i] : 0u;
      descs[i].window_len   = win_len ? win_len[i] : 0;
      descs[i].aux_len      = aux_len ? aux_len[i] : 0;
      descs[i].data_offset  = data_off ? data_off[i] : 0ull;
      descs[i].param_offset = param_off ? param_off[i] : 0ull;
      descs[i].payload_hint = payload_hint ? payload_hint[i] : 0u;
      descs[i].reserved     = 0u;
    }
    if (!WriteExact(h, descs.data(), static_cast<DWORD>(descs.size() * sizeof(descs[0])))) { ::CloseHandle(h); return -5; }
  }

  // Data blob (double[])
  if (dsize > 0 && data_blob) {
    if (!WriteExact(h, data_blob, static_cast<DWORD>(dsize))) { ::CloseHandle(h); return -5; }
  }
  // Params blob (uchar[])
  if (psize > 0 && params_blob) {
    if (!WriteExact(h, params_blob, static_cast<DWORD>(psize))) { ::CloseHandle(h); return -5; }
  }
  ::FlushFileBuffers(h);

  // Read BatchResponse (expect PENDING)
  wave_pipe::BatchResponse resp{};
  if (!ReadExact(h, &resp, sizeof(resp))) { ::CloseHandle(h); return -5; }
  ::CloseHandle(h);
  if (resp.magic != MESSAGE_MAGIC || resp.version != PROTOCOL_VERSION || resp.user_tag != tag) return -1;
  if (out_handle) *out_handle = tag;
  if (resp.status == static_cast<std::int32_t>(Status::PENDING)) return 0;
  return resp.status; // se não for PENDING, retorna erro/OK
}

__declspec(dllexport)
int AlglibBridge_BatchFetch(unsigned long long handle,
                            unsigned char* out_payload, int out_cap, int* out_len,
                            int* out_result_count,
                            int timeout_ms) {
  using namespace bridge;
  using namespace wave_pipe;
  HANDLE h = OpenPipe(timeout_ms);
  if (h == INVALID_HANDLE_VALUE) return -5;
  auto start = std::chrono::steady_clock::now();
  while (true) {
    FetchRequest fr{}; fr.magic = MESSAGE_MAGIC; fr.version = PROTOCOL_VERSION; fr.cmd = static_cast<std::uint32_t>(Command::BATCH_FETCH); fr.user_tag = handle;
    if (!WriteExact(h, &fr, sizeof(fr))) { ::CloseHandle(h); return -5; }
    ::FlushFileBuffers(h);

    BatchResponse resp{};
    if (!ReadExact(h, &resp, sizeof(resp))) { ::CloseHandle(h); return -5; }
    if (resp.magic != MESSAGE_MAGIC || resp.version != PROTOCOL_VERSION || resp.user_tag != handle) { ::CloseHandle(h); return -1; }
    if (resp.status == static_cast<std::int32_t>(Status::PENDING)) {
      if (timeout_ms >= 0) {
        auto now = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        if (ms > timeout_ms) { ::CloseHandle(h); return -3; }
      }
      ::Sleep(2);
      continue;
    }
    if (out_result_count) *out_result_count = static_cast<int>(resp.result_count);
    int need = (resp.payload_size > 0 && resp.payload_size < 0x7fffffff) ? static_cast<int>(resp.payload_size) : 0;
    if (out_len) *out_len = need;
    if (need > 0 && out_payload && out_cap > 0) {
      const int to_copy = (need < out_cap) ? need : out_cap;
      if (!ReadExact(h, out_payload, static_cast<DWORD>(to_copy))) { ::CloseHandle(h); return -5; }
      if (to_copy < need) {
        // drena o restante
        std::vector<unsigned char> tmp(need - to_copy);
        if (!ReadExact(h, tmp.data(), static_cast<DWORD>(tmp.size()))) { ::CloseHandle(h); return -5; }
      }
    } else if (need > 0) {
      // ler e descartar
      std::vector<unsigned char> tmp(need);
      if (!ReadExact(h, tmp.data(), static_cast<DWORD>(need))) { ::CloseHandle(h); return -5; }
    }
    ::CloseHandle(h);
    return resp.status; // 0 = OK; <0 = erro
  }
}

} // extern "C"
