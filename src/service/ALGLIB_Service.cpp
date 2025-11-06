#include "alglib_service_runtime.h"
#include "Logging.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#undef ERROR
#undef OK
#undef INVALID
#undef PENDING
#endif

namespace
{
// Versão hardcoded (não automática)
static const wchar_t* kServiceVersion = L"alglib_service 2025.11.06 00:00 BRT";
using gpu_service::FftResultPayload;
using gpu_service::ProcessResultPayload;
using gpu_service::BatchResultPayload;

struct ServiceConfig
  {
   std::wstring              pipe_name = L"\\\\.\\pipe\\alglib-wave_pipe";
   alglib_gpu::GpuConfig     gpu_config{alglib_gpu::BackendType::Auto, 0, 8, 1024, 5000};
   int                       worker_count = std::max(2u, std::thread::hardware_concurrency());
   alglib::logging::Level    log_level   = alglib::logging::Level::Info;
   std::wstring              log_path;
  };

std::wstring backend_to_label(alglib_gpu::BackendType type)
  {
   switch(type)
     {
      case alglib_gpu::BackendType::Cuda:      return L"CUDA";
      case alglib_gpu::BackendType::OpenCL:    return L"OpenCL";
      case alglib_gpu::BackendType::CpuFallback:return L"CPU";
      default:                                 return L"Desconhecido";
     }
  }

std::string narrow(const std::wstring& value)
  {
   if(value.empty())
      return {};
   const int size = WideCharToMultiByte(CP_UTF8, 0, value.c_str(), -1, nullptr, 0, nullptr, nullptr);
   if(size <= 0)
      return {};
   std::string result(static_cast<std::size_t>(size - 1), '\0');
   WideCharToMultiByte(CP_UTF8, 0, value.c_str(), -1, result.data(), size, nullptr, nullptr);
   return result;
  }

std::string narrow(const wchar_t* value)
  {
   if(value == nullptr)
      return {};
   return narrow(std::wstring(value));
  }

std::string get_host_name()
  {
   wchar_t buffer[MAX_COMPUTERNAME_LENGTH + 1];
   DWORD size = MAX_COMPUTERNAME_LENGTH + 1;
   if(!GetComputerNameW(buffer, &size))
      return {};
   return narrow(std::wstring(buffer, size));
  }

std::wstring widen(const std::string& value)
  {
#ifdef _WIN32
   if(value.empty())
      return std::wstring();
   const int len = MultiByteToWideChar(CP_UTF8, 0, value.c_str(), -1, nullptr, 0);
   if(len <= 0)
      return std::wstring();
   std::wstring result(static_cast<std::size_t>(len - 1), L'\0');
   MultiByteToWideChar(CP_UTF8, 0, value.c_str(), -1, result.data(), len);
   return result;
#else
   return std::wstring(value.begin(), value.end());
#endif
  }

const wchar_t* operation_to_string(wave_pipe::Operation op)
  {
   using wave_pipe::Operation;
   switch(op)
     {
      case Operation::PING:                    return L"PING";
      case Operation::FFT_COMPLEX_FORWARD:     return L"FFT_COMPLEX_FORWARD";
      case Operation::FFT_COMPLEX_INVERSE:     return L"FFT_COMPLEX_INVERSE";
      case Operation::FFT_REAL_FORWARD:        return L"FFT_REAL_FORWARD";
      case Operation::FFT_REAL_INVERSE:        return L"FFT_REAL_INVERSE";
      case Operation::FFT_SINE_TRANSFORM:      return L"FFT_SINE_TRANSFORM";
      case Operation::FFT_COSINE_TRANSFORM:    return L"FFT_COSINE_TRANSFORM";
      case Operation::FFT_TWO_REAL:            return L"FFT_TWO_REAL";
      case Operation::SPECTRAL_RESAMPLE:       return L"SPECTRAL_RESAMPLE";
      case Operation::SPECTRAL_ZERO_PAD:       return L"SPECTRAL_ZERO_PAD";
      case Operation::SPECTRAL_REMOVE_DC:      return L"SPECTRAL_REMOVE_DC";
      case Operation::SPECTRAL_FILTER_GAUSSIAN:return L"SPECTRAL_FILTER_GAUSSIAN";
      case Operation::SPECTRAL_FILTER_NOTCH:   return L"SPECTRAL_FILTER_NOTCH";
      case Operation::SPECTRAL_APPLY_MASK:     return L"SPECTRAL_APPLY_MASK";
      case Operation::SPECTRAL_DENOISE:        return L"SPECTRAL_DENOISE";
      case Operation::SPECTRAL_UPSCALE:        return L"SPECTRAL_UPSCALE";
      case Operation::SPECTRAL_DOWNSCALE:      return L"SPECTRAL_DOWNSCALE";
      case Operation::SPECTRAL_CONVOLUTION:    return L"SPECTRAL_CONVOLUTION";
      case Operation::SPECTRAL_CORRELATION:    return L"SPECTRAL_CORRELATION";
      case Operation::SPECTRAL_ANALYZE:        return L"SPECTRAL_ANALYZE";
      case Operation::SPECTRAL_DETECT_TRANSITIONS:return L"SPECTRAL_DETECT_TRANSITIONS";
      case Operation::SPECTRAL_INSTANT_METRICS:return L"SPECTRAL_INSTANT_METRICS";
      case Operation::SPECTRAL_PHASE_UNWRAP:   return L"SPECTRAL_PHASE_UNWRAP";
      case Operation::SSA_CREATE:              return L"SSA_CREATE";
      case Operation::SSA_DESTROY:             return L"SSA_DESTROY";
      case Operation::SSA_CONFIGURE:           return L"SSA_CONFIGURE";
      case Operation::SSA_CLEAR:               return L"SSA_CLEAR";
      case Operation::SSA_ADD_SEQUENCE:        return L"SSA_ADD_SEQUENCE";
      case Operation::SSA_APPEND_POINT:        return L"SSA_APPEND_POINT";
      case Operation::SSA_ANALYZE_LAST:        return L"SSA_ANALYZE_LAST";
      case Operation::SSA_FORECAST_LAST:       return L"SSA_FORECAST_LAST";
      default:                                 return L"DESCONHECIDO";
     }
  }

bool read_exact(HANDLE pipe, void* buffer, size_t bytes)
  {
   std::uint8_t* ptr = static_cast<std::uint8_t*>(buffer);
   size_t        remaining = bytes;
   while(remaining > 0)
     {
      DWORD read = 0;
      if(!ReadFile(pipe, ptr, static_cast<DWORD>(remaining), &read, nullptr))
         return false;
      if(read == 0)
         return false;
      remaining -= read;
      ptr += read;
     }
   return true;
  }

bool write_exact(HANDLE pipe, const void* buffer, size_t bytes)
  {
   const std::uint8_t* ptr = static_cast<const std::uint8_t*>(buffer);
   size_t remaining = bytes;
   while(remaining > 0)
     {
      DWORD written = 0;
      if(!WriteFile(pipe, ptr, static_cast<DWORD>(remaining), &written, nullptr))
         return false;
      if(written == 0)
         return false;
      remaining -= written;
      ptr += written;
     }
   return true;
  }

alglib_gpu::BackendType parse_backend(const std::wstring& value)
  {
   if(value == L"cuda")      return alglib_gpu::BackendType::Cuda;
   if(value == L"opencl")    return alglib_gpu::BackendType::OpenCL;
   if(value == L"cpu")       return alglib_gpu::BackendType::OpenCL; // mantém execução em GPU
   return alglib_gpu::BackendType::Auto;
  }

alglib::logging::Level parse_log_level(const std::wstring& value)
  {
   if(value == L"trace") return alglib::logging::Level::Trace;
   if(value == L"debug") return alglib::logging::Level::Debug;
   if(value == L"info")  return alglib::logging::Level::Info;
   if(value == L"warn")  return alglib::logging::Level::Warn;
   if(value == L"error") return alglib::logging::Level::Error;
   if(value == L"fatal") return alglib::logging::Level::Fatal;
   return alglib::logging::Level::Info;
  }

bool load_config(const std::wstring& path, ServiceConfig& cfg)
  {
   std::wifstream file(path);
   if(!file.is_open())
      return false;

   std::wstring line;
   while(std::getline(file, line))
     {
      if(line.empty() || line[0] == L'#' || line[0] == L';')
         continue;

      const auto eq = line.find(L'=');
      if(eq == std::wstring::npos)
         continue;

      const std::wstring key   = line.substr(0, eq);
      const std::wstring value = line.substr(eq + 1);

      if(key == L"PipeName")
         cfg.pipe_name = value;
      else if(key == L"Backend")
         cfg.gpu_config.backend = parse_backend(value);
      else if(key == L"DeviceIndex")
         cfg.gpu_config.device_index = std::stoi(value);
      else if(key == L"Streams")
         cfg.gpu_config.stream_count = std::max(1, std::stoi(value));
      else if(key == L"MinGpuWindow")
         cfg.gpu_config.min_gpu_window = std::max(64, std::stoi(value));
      else if(key == L"JobTimeoutMs")
         cfg.gpu_config.job_timeout_ms = std::max(1000, std::stoi(value));
      else if(key == L"Workers")
         cfg.worker_count = std::max(1, std::stoi(value));
      else if(key == L"LogLevel")
         cfg.log_level = parse_log_level(value);
      else if(key == L"LogFile")
         cfg.log_path = value;
     }
   return true;
  }

// FFT legacy helpers removidos

void send_process_pending(HANDLE pipe, std::uint64_t tag)
  {
   wave_pipe::ProcessResponse response{};
   response.magic            = wave_pipe::MESSAGE_MAGIC;
   response.version          = wave_pipe::PROTOCOL_VERSION;
   response.cmd              = static_cast<std::uint32_t>(wave_pipe::Command::PROCESS_FETCH);
   response.user_tag         = tag;
   response.status           = static_cast<std::int32_t>(wave_pipe::Status::PENDING);
   response.primary_count    = 0;
   response.secondary_count  = 0;
   response.cycle_count      = 0;
   response.instant_count    = 0;
   response.transition_count = 0;
   response.payload_size     = 0;
   write_exact(pipe, &response, sizeof(response));
  }

void send_process_response(HANDLE pipe, std::uint64_t tag, const ProcessResultPayload& payload)
  {
   wave_pipe::ProcessResponse header{};
   header.magic            = wave_pipe::MESSAGE_MAGIC;
   header.version          = wave_pipe::PROTOCOL_VERSION;
   header.cmd              = static_cast<std::uint32_t>(wave_pipe::Command::PROCESS_FETCH);
   header.user_tag         = tag;
   header.status           = static_cast<std::int32_t>(payload.status);
   header.primary_count    = static_cast<std::uint32_t>(payload.primary.size());
   header.secondary_count  = static_cast<std::uint32_t>(payload.secondary.size());
   header.cycle_count      = static_cast<std::uint32_t>(payload.cycles.size());
   header.instant_count    = static_cast<std::uint32_t>(payload.instant.size());
   header.transition_count = static_cast<std::uint32_t>(payload.transitions.size());
  header.payload_size     = static_cast<std::uint32_t>(payload.primary.size() * sizeof(double)
                                  + payload.secondary.size() * sizeof(double)
                                  + payload.cycles.size() * sizeof(wave_pipe::FftCyclePayload)
                                  + payload.instant.size() * sizeof(wave_pipe::InstantMetricsPayload)
                                  + payload.transitions.size() * sizeof(wave_pipe::TransitionPayload)
                                  + payload.metadata.size());

   write_exact(pipe, &header, sizeof(header));

   if(payload.status != wave_pipe::Status::OK)
      return;

   if(!payload.primary.empty())
      write_exact(pipe, payload.primary.data(), payload.primary.size() * sizeof(double));
   if(!payload.secondary.empty())
      write_exact(pipe, payload.secondary.data(), payload.secondary.size() * sizeof(double));
   if(!payload.cycles.empty())
      write_exact(pipe, payload.cycles.data(), payload.cycles.size() * sizeof(wave_pipe::FftCyclePayload));
  if(!payload.instant.empty())
     write_exact(pipe, payload.instant.data(), payload.instant.size() * sizeof(wave_pipe::InstantMetricsPayload));
  if(!payload.transitions.empty())
     write_exact(pipe, payload.transitions.data(), payload.transitions.size() * sizeof(wave_pipe::TransitionPayload));
  if(!payload.metadata.empty())
     write_exact(pipe, payload.metadata.data(), payload.metadata.size());
 }

void send_batch_pending(HANDLE pipe, std::uint64_t tag)
  {
   wave_pipe::BatchResponse response{};
   response.magic        = wave_pipe::MESSAGE_MAGIC;
   response.version      = wave_pipe::PROTOCOL_VERSION;
   response.cmd          = static_cast<std::uint32_t>(wave_pipe::Command::BATCH_FETCH);
   response.user_tag     = tag;
   response.status       = static_cast<std::int32_t>(wave_pipe::Status::PENDING);
   response.result_count = 0;
   response.payload_size = 0;
   response.reserved     = 0;
   write_exact(pipe, &response, sizeof(response));
  }

void send_batch_response(HANDLE pipe, std::uint64_t tag, const BatchResultPayload& payload)
  {
   wave_pipe::BatchResponse header{};
   header.magic        = wave_pipe::MESSAGE_MAGIC;
   header.version      = wave_pipe::PROTOCOL_VERSION;
   header.cmd          = static_cast<std::uint32_t>(wave_pipe::Command::BATCH_FETCH);
   header.user_tag     = tag;
   header.status       = static_cast<std::int32_t>(payload.status);
   header.result_count = static_cast<std::uint32_t>(payload.buffer.empty() ? 0 : 1);
   header.payload_size = static_cast<std::uint64_t>(payload.buffer.size());
   header.reserved     = 0;

   write_exact(pipe, &header, sizeof(header));
   if(payload.status != wave_pipe::Status::OK)
      return;
   if(!payload.buffer.empty())
      write_exact(pipe, payload.buffer.data(), payload.buffer.size());
  }
} // namespace

#ifdef _WIN32
int wmain(int argc, wchar_t* argv[])
  {
   ServiceConfig config;
   std::wstring  config_path;

   for(int i = 1; i < argc; ++i)
     {
      std::wstring arg = argv[i];
      if(arg == L"--config" && (i + 1) < argc)
         config_path = argv[++i];
      else if(arg == L"--pipe" && (i + 1) < argc)
         config.pipe_name = argv[++i];
      else if(arg == L"--device" && (i + 1) < argc)
         config.gpu_config.device_index = std::stoi(argv[++i]);
      else if(arg == L"--backend" && (i + 1) < argc)
         config.gpu_config.backend = parse_backend(argv[++i]);
      else if(arg == L"--workers" && (i + 1) < argc)
         config.worker_count = std::max(1, std::stoi(argv[++i]));
      else if(arg == L"--log-level" && (i + 1) < argc)
         config.log_level = parse_log_level(argv[++i]);
      else if(arg == L"--log-file" && (i + 1) < argc)
         config.log_path = argv[++i];
     }

  if(!config_path.empty())
    {
     if(!load_config(config_path, config))
        ALGLIB_LOG_WARN(L"Não foi possível carregar arquivo de configuração: " << config_path);
    }

  if(!config.log_path.empty())
     alglib::logging::SetLogFile(config.log_path);

  alglib::logging::SetLogLevel(config.log_level);
  ALGLIB_LOG_INFO(L"Versão: " << kServiceVersion);

  const std::wstring active_log_file = alglib::logging::GetLogFile();
  ALGLIB_LOG_INFO(L"Arquivo de log: " << active_log_file);

  gpu_service::ServiceRuntime runtime;
  if(!runtime.start(config.gpu_config, config.worker_count))
     {
      ALGLIB_LOG_FATAL(L"Falha ao inicializar backend GPU.");
      return 1;
     }

  const auto active_backend = runtime.resolved_backend();
  ALGLIB_LOG_INFO(L"Backend ativo: " << backend_to_label(active_backend));
  ALGLIB_LOG_INFO(L"Pipe: " << config.pipe_name);
  ALGLIB_LOG_INFO(L"Workers: " << config.worker_count);

  const std::wstring host_label = widen(get_host_name());
  std::wcout << L"==============================================================\n";
  std::wcout << L"  ALGLIB Service" << std::endl;
  std::wcout << L"  Versão    : " << kServiceVersion << std::endl;
  std::wcout << L"  Backend   : " << backend_to_label(active_backend) << std::endl;
  std::wcout << L"  Pipe      : " << config.pipe_name << std::endl;
  std::wcout << L"  Workers   : " << config.worker_count << std::endl;
  if(!host_label.empty())
     std::wcout << L"  Host      : " << host_label << std::endl;
  std::wcout << L"  PID       : " << GetCurrentProcessId() << std::endl;
  std::wcout << L"==============================================================\n";
  std::wcout << L"Pressione Ctrl+C para encerrar." << std::endl;
  std::wcout << std::endl;

  gpu_service::ServiceRuntime::Metadata metadata;
  metadata.version = narrow(kServiceVersion);
  metadata.backend = narrow(backend_to_label(active_backend));
  metadata.driver  = metadata.backend;
  metadata.build_timestamp = std::string(__DATE__ " " __TIME__);
  metadata.host    = get_host_name();
  metadata.workers = config.worker_count;
  metadata.pid     = static_cast<std::uint64_t>(GetCurrentProcessId());
  runtime.set_metadata(metadata);

  bool running = true;
  std::wcout << L"[service] Aguardando conexões no pipe " << config.pipe_name << L"..." << std::endl;
  while(running)
    {
      // Segurança do pipe: DACL permissiva para permitir conexão de processos não elevados.
      SECURITY_ATTRIBUTES sa{};
      SECURITY_DESCRIPTOR sd{};
      InitializeSecurityDescriptor(&sd, SECURITY_DESCRIPTOR_REVISION);
      // NULL DACL (permissivo) – adequado para ambiente de desenvolvimento/local.
      SetSecurityDescriptorDacl(&sd, TRUE, nullptr, FALSE);
      sa.nLength = sizeof(sa);
      sa.bInheritHandle = FALSE;
      sa.lpSecurityDescriptor = &sd;

      HANDLE pipe = CreateNamedPipeW(config.pipe_name.c_str(),
                                     PIPE_ACCESS_DUPLEX,
                                     PIPE_TYPE_MESSAGE | PIPE_READMODE_MESSAGE | PIPE_WAIT,
                                     PIPE_UNLIMITED_INSTANCES,
                                     1 << 20,
                                     1 << 20,
                                     5000,
                                     &sa);
      if(pipe == INVALID_HANDLE_VALUE)
        {
         ALGLIB_LOG_ERROR(L"CreateNamedPipeW falhou para pipe " << config.pipe_name);
         break;
        }

      BOOL connected = ConnectNamedPipe(pipe, nullptr) ? TRUE : (GetLastError() == ERROR_PIPE_CONNECTED);
      if(!connected)
        {
         ALGLIB_LOG_WARN(L"Cliente não conseguiu conectar-se ao pipe");
         std::wcout << L"[service] Falha ao conectar cliente." << std::endl;
         CloseHandle(pipe);
         continue;
        }

      ALGLIB_LOG_INFO(L"Cliente conectado ao pipe " << config.pipe_name);
      std::wcout << L"[service] Cliente conectado." << std::endl;

      bool client_active = true;
      while(client_active)
        {
         wave_pipe::FetchRequest header{};
         if(!read_exact(pipe, &header, sizeof(header)))
           {
            ALGLIB_LOG_WARN(L"Leitura de cabeçalho falhou; encerrando cliente");
            client_active = false;
            continue;
           }

         if(header.magic != wave_pipe::MESSAGE_MAGIC || header.version != wave_pipe::PROTOCOL_VERSION)
           {
            ALGLIB_LOG_WARN(L"Cabeçalho inválido recebido; cliente será desconectado");
            client_active = false;
            continue;
           }

         const auto cmd = static_cast<wave_pipe::Command>(header.cmd);

         switch(cmd)
           {

            case wave_pipe::Command::PROCESS_SUBMIT:
              {
               wave_pipe::ProcessRequest req{};
               req.magic    = header.magic;
               req.version  = header.version;
               req.cmd      = header.cmd;
               req.user_tag = header.user_tag;
               if(!read_exact(pipe,
                              reinterpret_cast<std::uint8_t*>(&req) + sizeof(wave_pipe::FetchRequest),
                              sizeof(wave_pipe::ProcessRequest) - sizeof(wave_pipe::FetchRequest)))
                 {
                  client_active = false;
                  break;
                 }

               std::vector<double> primary;
               std::vector<double> secondary;
               std::vector<std::uint8_t> params;
               if(req.window_len > 0)
                 {
                  primary.resize(static_cast<std::size_t>(req.window_len));
                  if(!read_exact(pipe, primary.data(), primary.size() * sizeof(double)))
                    {
                     client_active = false;
                     break;
                    }
                 }
               if(req.aux_len > 0)
                 {
                  secondary.resize(static_cast<std::size_t>(req.aux_len));
                  if(!read_exact(pipe, secondary.data(), secondary.size() * sizeof(double)))
                    {
                     client_active = false;
                     break;
                    }
                 }
               if(req.param_size > 0)
                 {
                  params.resize(req.param_size);
                  if(!read_exact(pipe, params.data(), params.size()))
                    {
                     client_active = false;
                     break;
                    }
                 }

               ALGLIB_LOG_DEBUG(L"PROCESS_SUBMIT recebido tag=" << req.user_tag << L" window=" << req.window_len << L" aux=" << req.aux_len << L" params=" << req.param_size);
               const auto op_name = operation_to_string(static_cast<wave_pipe::Operation>(req.operation));
               std::wcout << L"[service] PROCESS_SUBMIT tag=" << req.user_tag << L" op=" << op_name
                          << L" window=" << req.window_len << L" aux=" << req.aux_len << std::endl;
               runtime.submit_process(req,
                                      primary.empty() ? nullptr : primary.data(),
                                      secondary.empty() ? nullptr : secondary.data(),
                                      params.empty() ? nullptr : params.data());

               send_process_pending(pipe, req.user_tag);
              }
              break;

            case wave_pipe::Command::PROCESS_FETCH:
              {
               ProcessResultPayload payload;
               auto status = runtime.fetch_process(header.user_tag, payload);
               if(status == wave_pipe::Status::PENDING)
                 send_process_pending(pipe, header.user_tag);
               else
                 send_process_response(pipe, header.user_tag, payload);
              }
              break;

            case wave_pipe::Command::BATCH_SUBMIT:
              {
               wave_pipe::BatchHeader batch{};
               batch.magic    = header.magic;
               batch.version  = header.version;
               batch.cmd      = header.cmd;
               batch.user_tag = header.user_tag;
               if(!read_exact(pipe,
                              reinterpret_cast<std::uint8_t*>(&batch) + sizeof(wave_pipe::FetchRequest),
                              sizeof(wave_pipe::BatchHeader) - sizeof(wave_pipe::FetchRequest)))
                 {
                  client_active = false;
                  break;
                 }

               std::vector<wave_pipe::BatchTaskDescriptor> tasks;
               if(batch.task_count > 0)
                 {
                  tasks.resize(batch.task_count);
                  if(!read_exact(pipe, tasks.data(), tasks.size() * sizeof(wave_pipe::BatchTaskDescriptor)))
                    {
                     client_active = false;
                     break;
                    }
                 }

               std::vector<std::uint8_t> data_blob;
               if(batch.data_size > 0)
                 {
                  const std::size_t bytes = static_cast<std::size_t>(batch.data_size);
                  data_blob.resize(bytes);
                  if(!read_exact(pipe, data_blob.data(), bytes))
                    {
                     client_active = false;
                     break;
                    }
                 }

               std::vector<std::uint8_t> params_blob;
               if(batch.param_size > 0)
                 {
                  const std::size_t bytes = static_cast<std::size_t>(batch.param_size);
                  params_blob.resize(bytes);
                  if(!read_exact(pipe, params_blob.data(), bytes))
                    {
                     client_active = false;
                     break;
                    }
                 }

               ALGLIB_LOG_DEBUG(L"BATCH_SUBMIT recebido tag=" << batch.user_tag << L" tasks=" << batch.task_count);
               runtime.submit_batch(batch,
                                     tasks.empty() ? nullptr : tasks.data(),
                                     data_blob.empty() ? nullptr : data_blob.data(),
                                     params_blob.empty() ? nullptr : params_blob.data());

               send_batch_pending(pipe, batch.user_tag);
              }
              break;

            case wave_pipe::Command::BATCH_FETCH:
              {
               BatchResultPayload payload;
               auto status = runtime.fetch_batch(header.user_tag, payload);
               if(status == wave_pipe::Status::PENDING)
                 send_batch_pending(pipe, header.user_tag);
               else
                 send_batch_response(pipe, header.user_tag, payload);
              }
              break;

            default:
              send_pending_ack(pipe, header);
              break;
           }
        }

      DisconnectNamedPipe(pipe);
      CloseHandle(pipe);
      ALGLIB_LOG_INFO(L"Cliente desconectado.");
      std::wcout << L"[service] Cliente desconectado." << std::endl;
      std::wcout << L"[service] Aguardando conexões..." << std::endl;
     }

   runtime.stop();
   return 0;
  }
#else
int main()
  {
   std::cerr << "alglib_service só é suportado em Windows." << std::endl;
   return 1;
  }
#endif
