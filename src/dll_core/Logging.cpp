#include "Logging.h"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#ifndef _WIN32
#include <unistd.h>
#endif

namespace alglib
{
namespace logging
{

namespace
{
Level g_level = Level::Info;
std::wstring g_log_file;
std::mutex g_mutex;
bool g_initialized = false;

std::wstring WideFromUtf8(const std::string& text)
  {
#ifdef _WIN32
   if(text.empty())
      return std::wstring();
   const int count = MultiByteToWideChar(CP_UTF8, 0, text.c_str(), -1, nullptr, 0);
   if(count <= 0)
      return std::wstring();
   std::wstring result(static_cast<std::size_t>(count - 1), L'\0');
   MultiByteToWideChar(CP_UTF8, 0, text.c_str(), -1, result.data(), count);
   return result;
#else
   return std::wstring(text.begin(), text.end());
#endif
  }

std::wstring DefaultLogFile()
  {
#ifdef _WIN32
   wchar_t buffer[MAX_PATH];
   const DWORD len = GetModuleFileNameW(nullptr, buffer, MAX_PATH);
   std::filesystem::path path;
   if(len > 0)
      path = std::filesystem::path(buffer).parent_path();
   else
      path = std::filesystem::temp_directory_path();
   path /= L"alglib_service.log";
   return path.wstring();
#else
   return L"alglib_service.log";
#endif
  }

std::string ToUtf8(const std::wstring& text)
  {
#ifdef _WIN32
   if(text.empty())
      return std::string();
   const int count = WideCharToMultiByte(CP_UTF8, 0, text.c_str(), (int)text.size(), nullptr, 0, nullptr, nullptr);
   if(count <= 0)
      return std::string();
   std::string result(static_cast<std::size_t>(count), '\0');
   WideCharToMultiByte(CP_UTF8, 0, text.c_str(), (int)text.size(), result.data(), count, nullptr, nullptr);
   return result;
#else
   return std::string(text.begin(), text.end());
#endif
  }

std::string LevelToString(Level level)
  {
   switch(level)
     {
      case Level::Trace: return "TRACE";
      case Level::Debug: return "DEBUG";
      case Level::Info:  return "INFO";
      case Level::Warn:  return "WARN";
      case Level::Error: return "ERROR";
      case Level::Fatal: return "FATAL";
      default:           return "UNKNOWN";
     }
  }

void EnsureInitialized()
  {
   if(g_initialized)
      return;

   g_log_file = DefaultLogFile();
   g_initialized = true;
  }

std::wstring FromNarrow(const char* text)
  {
   if(text == nullptr)
      return std::wstring();
   return WideFromUtf8(text);
  }

} // namespace

void SetLogLevel(Level level)
  {
   g_level = level;
  }

Level GetLogLevel()
  {
   return g_level;
  }

void SetLogFile(const std::wstring& path)
  {
   std::lock_guard<std::mutex> lock(g_mutex);
   g_log_file = path;
   g_initialized = true;
  }

std::wstring GetLogFile()
  {
   std::lock_guard<std::mutex> lock(g_mutex);
   EnsureInitialized();
   return g_log_file;
  }

void Log(Level level,
         const char* file,
         int line,
         const char* function,
         const std::wstring& message)
  {
   if(level < g_level)
      return;

   std::lock_guard<std::mutex> lock(g_mutex);
   EnsureInitialized();

   std::ofstream stream(ToUtf8(g_log_file), std::ios::app | std::ios::binary);
   if(!stream.is_open())
      return;

   const auto now = std::chrono::system_clock::now();
   const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
   const std::time_t tt = std::chrono::system_clock::to_time_t(now);
   std::tm tm_snapshot;
#ifdef _WIN32
   localtime_s(&tm_snapshot, &tt);
#else
   tm_snapshot = *std::localtime(&tt);
#endif

#ifdef _WIN32
   const DWORD pid = GetCurrentProcessId();
   const DWORD tid = GetCurrentThreadId();
#else
   const auto pid = static_cast<long>(::getpid());
   const auto tid = std::hash<std::thread::id>{}(std::this_thread::get_id());
#endif

   std::ostringstream oss;
   oss << std::put_time(&tm_snapshot, "%Y-%m-%d %H:%M:%S")
       << '.' << std::setw(3) << std::setfill('0') << ms.count()
       << ' ' << '[' << LevelToString(level) << ']'
       << " [PID:" << pid << "]"
       << " [TID:" << tid << "]"
       << " [" << file << ':' << line << ']'
       << ' ' << function << " - " << ToUtf8(message);

  stream << oss.str() << std::endl;
  std::cout << oss.str() << std::endl;
  std::cout.flush();
 }

void Log(Level level,
         const char* file,
         int line,
         const char* function,
         const std::string& message)
  {
   Log(level, file, line, function, WideFromUtf8(message));
  }

} // namespace logging
} // namespace alglib
