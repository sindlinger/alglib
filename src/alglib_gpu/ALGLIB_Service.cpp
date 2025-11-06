// alglib_service.cpp — serviço principal (pipe único) para Via A
// Atenção: versão hardcoded (exigência do projeto)

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>

#include "ap.h"
#include "alglib_pipe_messages.h"
#include "alglib_service_runtime.h"
#include "GpuBackend.h"
#include "Logging.h"

#define NOMINMAX
#include <windows.h>
#undef OK
#undef ERROR
#undef INVALID
#undef PENDING

namespace
{
constexpr wchar_t kPipeName[] = L"\\\\.\\pipe\\alglib-wave_pipe"; // sem "2"
constexpr DWORD   kPipeOutBuffer = 1 << 20;
constexpr DWORD   kPipeInBuffer  = 1 << 20;
constexpr int     kDefaultWorkers = 24;

// Versão hardcoded conforme solicitação (não automatizar)
// ATENÇÃO: Atualize manualmente sempre que fizer um build público.
constexpr const char* kVersionString = "alglib_service 2025.11.06 03:50 BRT";

bool ReadExact(HANDLE pipe, void* buffer, DWORD bytes)
{
    auto* dst = static_cast<std::uint8_t*>(buffer);
    DWORD total = 0;
    while (total < bytes)
    {
        DWORD chunk = 0;
        if (!ReadFile(pipe, dst + total, bytes - total, &chunk, nullptr) || chunk == 0)
            return false;
        total += chunk;
    }
    return true;
}

bool WriteExact(HANDLE pipe, const void* buffer, DWORD bytes)
{
    auto* src = static_cast<const std::uint8_t*>(buffer);
    DWORD total = 0;
    while (total < bytes)
    {
        DWORD chunk = 0;
        if (!WriteFile(pipe, src + total, bytes - total, &chunk, nullptr) || chunk == 0)
            return false;
        total += chunk;
    }
    return true;
}

void PrintHeader()
{
    std::printf("=============================================================\n");
    std::printf("  alglib_service\n");
    std::printf("  Versão    : %s\n", kVersionString);
    std::printf("  Backend   : (dynamic)\n");
    std::printf("  Pipe      : \\.\\pipe\\alglib-wave_pipe\n");
    std::printf("  Workers   : %d\n", kDefaultWorkers);
    std::printf("  Host      : (local)\n");
    std::printf("=============================================================\n");
    std::printf("Pressione Ctrl+C para encerrar.\n\n");
    std::fflush(stdout);
}

void ServeOnce(gpu_service::ServiceRuntime& runtime)
{
    ALGLIB_LOG_INFO(L"ServeOnce: criando named pipe em '" << kPipeName << L"'");
    HANDLE pipe = CreateNamedPipeW(kPipeName,
                                   PIPE_ACCESS_DUPLEX,
                                   PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
                                   1,
                                   kPipeOutBuffer,
                                   kPipeInBuffer,
                                   0,
                                   nullptr);
    if (pipe == INVALID_HANDLE_VALUE)
    {
        ALGLIB_LOG_ERROR(L"CreateNamedPipeW falhou: " << (int)GetLastError());
        Sleep(250);
        return;
    }

    BOOL connected = ConnectNamedPipe(pipe, nullptr) ? TRUE : (GetLastError() == ERROR_PIPE_CONNECTED);
    if (!connected)
    {
        ALGLIB_LOG_WARN(L"ConnectNamedPipe: cliente não conectado (erro=" << (int)GetLastError() << L")");
        CloseHandle(pipe);
        Sleep(50);
        return;
    }
    ALGLIB_LOG_INFO(L"Cliente conectado ao pipe");

    // Loop por conexão: responde PROCESS_SUBMIT/FETCH e BATCH_SUBMIT/FETCH
    while (true)
    {
        // Peek para detectar desconexão limpa
        DWORD avail = 0;
        if (!PeekNamedPipe(pipe, nullptr, 0, nullptr, &avail, nullptr))
        {
            ALGLIB_LOG_INFO(L"Cliente desconectou (PeekNamedPipe=false)");
            break;
        }

        // Escolha de cabeçalho: ProcessRequest (36 bytes + 12) = 36? Não — usamos struct completa
        // Lemos magic/version/cmd para decidir o tipo
        wave_pipe::FetchRequest peek{}; // usa layout compatível no início
        if (!ReadExact(pipe, &peek, sizeof(peek)))
        {
            ALGLIB_LOG_ERROR(L"Falha lendo FetchRequest/peek (ReadExact=false)");
            break;
        }

        if (peek.magic != wave_pipe::MESSAGE_MAGIC || peek.version != wave_pipe::PROTOCOL_VERSION)
        {
            ALGLIB_LOG_ERROR(L"Mensagem inválida (magic/version)");
            break;
        }

        const auto cmd = static_cast<wave_pipe::Command>(peek.cmd);
        if (cmd == wave_pipe::Command::PROCESS_FETCH)
        {
            ALGLIB_LOG_DEBUG(L"Recebido PROCESS_FETCH (tag=" << peek.user_tag << L")");
            // Já temos o FetchRequest completo em 'peek'
            gpu_service::ProcessResultPayload payload{};
            auto status = runtime.fetch_process(peek.user_tag, payload);

            wave_pipe::ProcessResponse resp{};
            resp.magic = wave_pipe::MESSAGE_MAGIC;
            resp.version = wave_pipe::PROTOCOL_VERSION;
            resp.cmd = static_cast<std::uint32_t>(wave_pipe::Command::PROCESS_FETCH);
            resp.user_tag = peek.user_tag;
            resp.status = static_cast<std::int32_t>(status);

            std::vector<std::uint8_t> bytes;
            if (status == wave_pipe::Status::OK)
            {
                resp.primary_count    = static_cast<std::uint32_t>(payload.primary.size());
                resp.secondary_count  = static_cast<std::uint32_t>(payload.secondary.size());
                resp.cycle_count      = static_cast<std::uint32_t>(payload.cycles.size());
                resp.instant_count    = static_cast<std::uint32_t>(payload.instant.size());
                resp.transition_count = static_cast<std::uint32_t>(payload.transitions.size());

                const std::size_t payload_bytes =
                    payload.primary.size()    * sizeof(double) +
                    payload.secondary.size()  * sizeof(double) +
                    payload.cycles.size()     * sizeof(wave_pipe::FftCyclePayload) +
                    payload.instant.size()    * sizeof(wave_pipe::InstantMetricsPayload) +
                    payload.transitions.size()* sizeof(wave_pipe::TransitionPayload) +
                    payload.metadata.size();

                resp.payload_size = static_cast<std::uint32_t>(payload_bytes);
                bytes.reserve(payload_bytes);

                auto append_raw = [&bytes](const void* data, std::size_t n){
                    const auto* p = static_cast<const std::uint8_t*>(data);
                    bytes.insert(bytes.end(), p, p + n);
                };
                if (!payload.primary.empty())     append_raw(payload.primary.data(),     payload.primary.size()     * sizeof(double));
                if (!payload.secondary.empty())   append_raw(payload.secondary.data(),   payload.secondary.size()   * sizeof(double));
                if (!payload.cycles.empty())      append_raw(payload.cycles.data(),      payload.cycles.size()      * sizeof(wave_pipe::FftCyclePayload));
                if (!payload.instant.empty())     append_raw(payload.instant.data(),     payload.instant.size()     * sizeof(wave_pipe::InstantMetricsPayload));
                if (!payload.transitions.empty()) append_raw(payload.transitions.data(), payload.transitions.size() * sizeof(wave_pipe::TransitionPayload));
                if (!payload.metadata.empty())    append_raw(payload.metadata.data(),    payload.metadata.size());
            }

            if (!WriteExact(pipe, &resp, sizeof(resp)))
            {
                ALGLIB_LOG_ERROR(L"Falha escrevendo ProcessResponse header (FETCH)");
                break;
            }
            if (resp.payload_size > 0)
            {
                if (!WriteExact(pipe, bytes.data(), static_cast<DWORD>(bytes.size())))
                {
                    ALGLIB_LOG_ERROR(L"Falha escrevendo payload (FETCH)");
                    break;
                }
            }
            FlushFileBuffers(pipe);
        }
        else if (cmd == wave_pipe::Command::PROCESS_SUBMIT)
        {
            ALGLIB_LOG_DEBUG(L"Recebido PROCESS_SUBMIT (tag=" << peek.user_tag << L")");
            // Reinterpretar o buffer como ProcessRequest (os primeiros campos já lidos são compatíveis)
            wave_pipe::ProcessRequest req{};
            req.magic   = peek.magic;
            req.version = peek.version;
            req.cmd     = peek.cmd;
            req.user_tag= peek.user_tag;

            // Ler o restante do ProcessRequest
            std::uint8_t rest[sizeof(wave_pipe::ProcessRequest) - sizeof(wave_pipe::FetchRequest)]{};
            if (!ReadExact(pipe, rest, sizeof(rest)))
            {
                ALGLIB_LOG_ERROR(L"Falha lendo restante do ProcessRequest");
                break;
            }
            std::memcpy(reinterpret_cast<std::uint8_t*>(&req) + sizeof(wave_pipe::FetchRequest), rest, sizeof(rest));

            std::vector<std::uint8_t> params(req.param_size);
            if (req.param_size > 0)
            {
                if (!ReadExact(pipe, params.data(), req.param_size))
                {
                    ALGLIB_LOG_ERROR(L"Falha lendo params (size=" << (int)req.param_size << L")");
                    break;
                }
            }
            std::vector<double> primary(static_cast<std::size_t>(std::max(0, req.window_len)));
            if (req.window_len > 0)
            {
                if (!ReadExact(pipe, primary.data(), static_cast<DWORD>(primary.size() * sizeof(double))))
                {
                    ALGLIB_LOG_ERROR(L"Falha lendo primary (len=" << req.window_len << L")");
                    break;
                }
            }
            std::vector<double> secondary(static_cast<std::size_t>(std::max(0, req.aux_len)));
            if (req.aux_len > 0)
            {
                if (!ReadExact(pipe, secondary.data(), static_cast<DWORD>(secondary.size() * sizeof(double))))
                {
                    ALGLIB_LOG_ERROR(L"Falha lendo secondary (len=" << req.aux_len << L")");
                    break;
                }
            }

            // Submete no runtime (assíncrono)
            ALGLIB_LOG_INFO(L"SUBMIT op=" << req.operation << L" win=" << req.window_len << L" aux=" << req.aux_len << L" params=" << (int)req.param_size);
            runtime.submit_process(req,
                                   primary.empty()   ? nullptr : primary.data(),
                                   secondary.empty() ? nullptr : secondary.data(),
                                   params.empty()    ? nullptr : params.data());

            // Responde imediatamente com PENDING (cliente fará FETCH)
            wave_pipe::ProcessResponse resp{};
            resp.magic = wave_pipe::MESSAGE_MAGIC;
            resp.version = wave_pipe::PROTOCOL_VERSION;
            resp.cmd = static_cast<std::uint32_t>(wave_pipe::Command::PROCESS_SUBMIT);
            resp.user_tag = req.user_tag;
            resp.status = static_cast<std::int32_t>(wave_pipe::Status::PENDING);
            resp.primary_count = 0;
            resp.secondary_count = 0;
            resp.cycle_count = 0;
            resp.instant_count = 0;
            resp.transition_count = 0;
            resp.payload_size = 0;
            if (!WriteExact(pipe, &resp, sizeof(resp)))
            {
                ALGLIB_LOG_ERROR(L"Falha escrevendo ProcessResponse (PENDING)");
                break;
            }
            FlushFileBuffers(pipe);
        }
        else if (cmd == wave_pipe::Command::BATCH_SUBMIT)
        {
            ALGLIB_LOG_DEBUG(L"Recebido BATCH_SUBMIT (tag=" << peek.user_tag << L")");
            // Ler cabeçalho completo
            wave_pipe::BatchHeader hdr{};
            hdr.magic = peek.magic; hdr.version = peek.version; hdr.cmd = peek.cmd; hdr.user_tag = peek.user_tag;
            std::uint8_t rest[sizeof(wave_pipe::BatchHeader) - sizeof(wave_pipe::FetchRequest)]{};
            if(!ReadExact(pipe, rest, sizeof(rest))) { ALGLIB_LOG_ERROR(L"Falha lendo restante do BatchHeader"); break; }
            std::memcpy(reinterpret_cast<std::uint8_t*>(&hdr) + sizeof(wave_pipe::FetchRequest), rest, sizeof(rest));

            std::vector<wave_pipe::BatchTaskDescriptor> tasks(hdr.task_count);
            if(hdr.task_count > 0)
            {
                if(!ReadExact(pipe, tasks.data(), static_cast<DWORD>(tasks.size() * sizeof(wave_pipe::BatchTaskDescriptor)))) { ALGLIB_LOG_ERROR(L"Falha lendo descriptors do batch"); break; }
            }

            std::vector<std::uint8_t> data_blob(static_cast<std::size_t>(hdr.data_size));
            if(hdr.data_size > 0)
            {
                if(!ReadExact(pipe, data_blob.data(), static_cast<DWORD>(data_blob.size()))) { ALGLIB_LOG_ERROR(L"Falha lendo data_blob do batch"); break; }
            }

            std::vector<std::uint8_t> params_blob(static_cast<std::size_t>(hdr.param_size));
            if(hdr.param_size > 0)
            {
                if(!ReadExact(pipe, params_blob.data(), static_cast<DWORD>(params_blob.size()))) { ALGLIB_LOG_ERROR(L"Falha lendo params_blob do batch"); break; }
            }

            runtime.submit_batch(hdr,
                                 tasks.empty() ? nullptr : tasks.data(),
                                 data_blob.empty() ? nullptr : data_blob.data(),
                                 params_blob.empty() ? nullptr : params_blob.data());

            wave_pipe::BatchResponse resp{};
            resp.magic = wave_pipe::MESSAGE_MAGIC;
            resp.version = wave_pipe::PROTOCOL_VERSION;
            resp.cmd = static_cast<std::uint32_t>(wave_pipe::Command::BATCH_SUBMIT);
            resp.user_tag = hdr.user_tag;
            resp.status = static_cast<std::int32_t>(wave_pipe::Status::PENDING);
            resp.result_count = 0;
            resp.payload_size = 0;
            resp.reserved = 0;
            if(!WriteExact(pipe, &resp, sizeof(resp))) { ALGLIB_LOG_ERROR(L"Falha escrevendo BatchResponse (PENDING)"); break; }
            FlushFileBuffers(pipe);
        }
        else if (cmd == wave_pipe::Command::BATCH_FETCH)
        {
            ALGLIB_LOG_DEBUG(L"Recebido BATCH_FETCH (tag=" << peek.user_tag << L")");
            gpu_service::BatchResultPayload payload{};
            auto status = runtime.fetch_batch(peek.user_tag, payload);

            wave_pipe::BatchResponse resp{};
            resp.magic = wave_pipe::MESSAGE_MAGIC;
            resp.version = wave_pipe::PROTOCOL_VERSION;
            resp.cmd = static_cast<std::uint32_t>(wave_pipe::Command::BATCH_FETCH);
            resp.user_tag = peek.user_tag;
            resp.status = static_cast<std::int32_t>(status);
            resp.result_count = (status == wave_pipe::Status::OK) ? 1u : 0u;
            resp.payload_size = (status == wave_pipe::Status::OK) ? static_cast<std::uint64_t>(payload.buffer.size()) : 0ull;
            resp.reserved = 0;

            if(!WriteExact(pipe, &resp, sizeof(resp))) { ALGLIB_LOG_ERROR(L"Falha escrevendo BatchResponse header (FETCH)"); break; }
            if(resp.payload_size > 0)
            {
                if(!WriteExact(pipe, payload.buffer.data(), static_cast<DWORD>(payload.buffer.size()))) { ALGLIB_LOG_ERROR(L"Falha escrevendo payload do batch"); break; }
            }
            FlushFileBuffers(pipe);
        }
        else
        {
            ALGLIB_LOG_WARN(L"Comando desconhecido recebido: " << (int)peek.cmd);
            break; // comando desconhecido
        }
    }

    FlushFileBuffers(pipe);
    DisconnectNamedPipe(pipe);
    CloseHandle(pipe);
    ALGLIB_LOG_INFO(L"Conexão encerrada");
}

} // namespace

int wmain()
{
    SetConsoleOutputCP(CP_UTF8);
    PrintHeader();
    // Inicializa logging para arquivo padrão ao lado do executável
    using alglib::logging::SetLogLevel; using alglib::logging::Level; using alglib::logging::GetLogFile;
    SetLogLevel(Level::Info);
    ALGLIB_LOG_INFO(L"alglib_service iniciado. Log em: " << GetLogFile());

    // Inicializa runtime/metadata
    gpu_service::ServiceRuntime runtime;
    alglib_gpu::GpuConfig cfg{};
    runtime.start(cfg, kDefaultWorkers);

    gpu_service::ServiceRuntime::Metadata meta;
    meta.version = kVersionString;
    meta.backend = "CUDA"; // O runtime preenche resolved_backend; aqui fica apenas informativo
    meta.driver  = "";
    meta.build_timestamp = __DATE__ " " __TIME__;
    meta.host = "";
    meta.workers = kDefaultWorkers;
    meta.pid = static_cast<std::uint64_t>(GetCurrentProcessId());
    runtime.set_metadata(meta);

    // Loop principal (uma conexão por vez; backend é multi-thread)
    while (true)
    {
        ServeOnce(runtime);
        Sleep(25);
    }

    return 0;
}
