#define NOMINMAX
#include <windows.h>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>

namespace
{
constexpr wchar_t kPipeName[] = L"\\\\.\\pipe\\alglib-wave_pipe";
constexpr uint32_t kMessageMagic = 0x4C574650; // 'P' 'F' 'W' 'L'
constexpr uint32_t kProtocolVersion = 1;
constexpr uint32_t kCmdProcessSubmit = 20;
constexpr uint32_t kCmdProcessFetch  = 21;
constexpr int32_t  kStatusOk          = 0;
constexpr int32_t  kStatusInvalid     = -2;
constexpr int32_t  kStatusIoError     = -5;

#pragma pack(push, 1)
struct RequestHeader
{
    uint32_t magic;
    uint32_t version;
    uint32_t cmd;
    uint64_t tag;
    uint32_t operation;
    uint32_t reserved;
    int32_t  primary_len;
    int32_t  secondary_len;
    uint32_t param_len;
};

struct ResponseHeader
{
    uint32_t magic;
    uint32_t version;
    uint32_t cmd;
    uint64_t tag;
    int32_t  status;
    uint32_t primary_count;
    uint32_t secondary_count;
    uint32_t cycles_count;
    uint32_t instant_count;
    uint32_t transition_count;
    uint32_t payload_size;
};
#pragma pack(pop)

bool ReadExact(HANDLE pipe, void *buffer, DWORD bytes)
{
    auto *dst = static_cast<uint8_t *>(buffer);
    DWORD total = 0;
    while (total < bytes)
    {
        DWORD chunk = 0;
        if (!ReadFile(pipe, dst + total, bytes - total, &chunk, nullptr))
        {
            return false;
        }
        if (chunk == 0)
        {
            return false;
        }
        total += chunk;
    }
    return true;
}

bool WriteExact(HANDLE pipe, const void *buffer, DWORD bytes)
{
    auto *src = static_cast<const uint8_t *>(buffer);
    DWORD total = 0;
    while (total < bytes)
    {
        DWORD chunk = 0;
        if (!WriteFile(pipe, src + total, bytes - total, &chunk, nullptr))
        {
            return false;
        }
        if (chunk == 0)
        {
            return false;
        }
        total += chunk;
    }
    return true;
}

void Log(const char *level, const std::string &message)
{
    SYSTEMTIME st{};
    GetLocalTime(&st);
    std::printf("[%04d-%02d-%02d %02d:%02d:%02d.%03d][DiagService][%s] %s\n",
                st.wYear,
                st.wMonth,
                st.wDay,
                st.wHour,
                st.wMinute,
                st.wSecond,
                st.wMilliseconds,
                level,
                message.c_str());
    std::fflush(stdout);
}

std::string ErrorMessage(DWORD err)
{
    LPWSTR buffer = nullptr;
    DWORD len = FormatMessageW(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                               nullptr,
                               err,
                               0,
                               reinterpret_cast<LPWSTR>(&buffer),
                               0,
                               nullptr);
    std::string result;
    if (len && buffer)
    {
        int utf8_len = WideCharToMultiByte(CP_UTF8, 0, buffer, -1, nullptr, 0, nullptr, nullptr);
        if (utf8_len > 0)
        {
            result.resize(static_cast<size_t>(utf8_len));
            WideCharToMultiByte(CP_UTF8, 0, buffer, -1, result.data(), utf8_len, nullptr, nullptr);
        }
        LocalFree(buffer);
    }
    if (result.empty())
    {
        result = "erro desconhecido";
    }
    return result;
}

void ProcessPipe()
{
    HANDLE pipe = CreateNamedPipeW(kPipeName,
                                   PIPE_ACCESS_DUPLEX,
                                   PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
                                   1,
                                   1 << 20,
                                   1 << 20,
                                   0,
                                   nullptr);

    if (pipe == INVALID_HANDLE_VALUE)
    {
        DWORD err = GetLastError();
        Log("ERROR", "CreateNamedPipe falhou: " + ErrorMessage(err));
        Sleep(1000);
        return;
    }

    Log("INFO", "Aguardando conexão no pipe alglib-wave_pipe...");

    BOOL connected = ConnectNamedPipe(pipe, nullptr) ? TRUE : (GetLastError() == ERROR_PIPE_CONNECTED);
    if (!connected)
    {
        DWORD err = GetLastError();
        Log("ERROR", "ConnectNamedPipe falhou: " + ErrorMessage(err));
        CloseHandle(pipe);
        Sleep(200);
        return;
    }

    Log("INFO", "Cliente conectado.");

    bool running = true;
    while (running)
    {
        RequestHeader req{};
        DWORD available = 0;
        if (!PeekNamedPipe(pipe, nullptr, 0, nullptr, &available, nullptr))
        {
            DWORD err = GetLastError();
            Log("WARN", "PeekNamedPipe falhou: " + ErrorMessage(err));
            break;
        }
        if (available == 0)
        {
            // Bloqueia até receber header.
        }

        if (!ReadExact(pipe, &req, sizeof(req)))
        {
            DWORD err = GetLastError();
            if (err == ERROR_BROKEN_PIPE || err == ERROR_NO_DATA)
            {
                Log("INFO", "Cliente desconectou.");
            }
            else
            {
                Log("ERROR", "ReadExact falhou no header: " + ErrorMessage(err));
            }
            break;
        }

        if (req.magic != kMessageMagic || req.version != kProtocolVersion)
        {
            Log("ERROR", "Mensagem inválida recebida (magic/version).");
            running = false;
            break;
        }

        if (req.cmd == kCmdProcessFetch)
        {
            Log("INFO", "Recebido FETCH (tag=" + std::to_string(req.tag) + "). Respondendo vazio.");
            ResponseHeader resp{};
            resp.magic = kMessageMagic;
            resp.version = kProtocolVersion;
            resp.cmd = kCmdProcessFetch;
            resp.tag = req.tag;
            resp.status = kStatusOk;
            resp.primary_count = 0;
            resp.secondary_count = 0;
            resp.cycles_count = 0;
            resp.instant_count = 0;
            resp.transition_count = 0;
            resp.payload_size = 0;
            WriteExact(pipe, &resp, sizeof(resp));
            FlushFileBuffers(pipe);
            continue;
        }

        if (req.cmd != kCmdProcessSubmit)
        {
            Log("WARN", "Comando desconhecido: " + std::to_string(req.cmd));
            running = false;
            break;
        }

        Log("INFO", "SUBMIT op=" + std::to_string(req.operation) +
                         ", primary=" + std::to_string(req.primary_len) +
                         ", secondary=" + std::to_string(req.secondary_len) +
                         ", params=" + std::to_string(req.param_len));

        int32_t primary_len = std::max(req.primary_len, 0);
        int32_t secondary_len = std::max(req.secondary_len, 0);

        std::vector<uint8_t> params(static_cast<size_t>(req.param_len));
        if (req.param_len > 0)
        {
            if (!ReadExact(pipe, params.data(), req.param_len))
            {
                Log("ERROR", "Falha lendo params.");
                break;
            }
        }

        std::vector<double> primary(static_cast<size_t>(primary_len));
        if (primary_len > 0)
        {
            if (!ReadExact(pipe, primary.data(), static_cast<DWORD>(primary.size() * sizeof(double))))
            {
                Log("ERROR", "Falha lendo primary payload.");
                break;
            }
        }

        std::vector<double> secondary(static_cast<size_t>(secondary_len));
        if (secondary_len > 0)
        {
            if (!ReadExact(pipe, secondary.data(), static_cast<DWORD>(secondary.size() * sizeof(double))))
            {
                Log("ERROR", "Falha lendo secondary payload.");
                break;
            }
        }

        // Para diagnóstico devolvemos os mesmos buffers.
        ResponseHeader resp{};
        resp.magic = kMessageMagic;
        resp.version = kProtocolVersion;
        resp.cmd = kCmdProcessSubmit;
        resp.tag = req.tag;
        resp.status = kStatusOk;
        resp.primary_count = static_cast<uint32_t>(primary.size());
        resp.secondary_count = static_cast<uint32_t>(secondary.size());
        resp.cycles_count = 0;
        resp.instant_count = 0;
        resp.transition_count = 0;
        resp.payload_size = static_cast<uint32_t>((primary.size() + secondary.size()) * sizeof(double));

        std::vector<double> out_payload;
        out_payload.reserve(primary.size() + secondary.size());
        out_payload.insert(out_payload.end(), primary.begin(), primary.end());
        out_payload.insert(out_payload.end(), secondary.begin(), secondary.end());

        if (!WriteExact(pipe, &resp, sizeof(resp)))
        {
            Log("ERROR", "Falha ao enviar header de resposta.");
            break;
        }
        if (!out_payload.empty())
        {
            if (!WriteExact(pipe, out_payload.data(), static_cast<DWORD>(out_payload.size() * sizeof(double))))
            {
                Log("ERROR", "Falha ao enviar payload de resposta.");
                break;
            }
        }
        FlushFileBuffers(pipe);
        Log("INFO", "Resposta enviada (tag=" + std::to_string(req.tag) + ").");
    }

    FlushFileBuffers(pipe);
    DisconnectNamedPipe(pipe);
    CloseHandle(pipe);
    Log("INFO", "Conexão encerrada. Reiniciando loop.");
}

} // namespace

int wmain()
{
    SetConsoleOutputCP(CP_UTF8);
    Log("INFO", "ALGLIB Diagnostic Service inicializado.");

    while (true)
    {
        ProcessPipe();
        Sleep(200);
    }

    return 0;
}
