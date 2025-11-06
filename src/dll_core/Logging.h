#pragma once

#include <string>

namespace alglib
{
namespace logging
{

enum class Level
  {
   Trace = 0,
   Debug,
   Info,
   Warn,
   Error,
   Fatal
  };

void SetLogLevel(Level level);
void SetLogFile(const std::wstring& path);
void Log(Level level,
         const char* file,
         int line,
         const char* function,
         const std::wstring& message);
void Log(Level level,
         const char* file,
         int line,
         const char* function,
         const std::string& message);

Level GetLogLevel();

std::wstring GetLogFile();

} // namespace logging
} // namespace alglib

#include <sstream>

#define ALGLIB_LOG_STREAM(level, stream_expr)                                        \
  do                                                                                \
    {                                                                               \
     std::wostringstream _alglib_log_stream;                                        \
     _alglib_log_stream << stream_expr;                                             \
     ::alglib::logging::Log(level, __FILE__, __LINE__, __FUNCTION__, _alglib_log_stream.str()); \
    }                                                                               \
  while(false)

#define ALGLIB_LOG_TRACE(message_expr) ALGLIB_LOG_STREAM(::alglib::logging::Level::Trace, message_expr)
#define ALGLIB_LOG_DEBUG(message_expr) ALGLIB_LOG_STREAM(::alglib::logging::Level::Debug, message_expr)
#define ALGLIB_LOG_INFO(message_expr)  ALGLIB_LOG_STREAM(::alglib::logging::Level::Info,  message_expr)
#define ALGLIB_LOG_WARN(message_expr)  ALGLIB_LOG_STREAM(::alglib::logging::Level::Warn,  message_expr)
#define ALGLIB_LOG_ERROR(message_expr) ALGLIB_LOG_STREAM(::alglib::logging::Level::Error, message_expr)
#define ALGLIB_LOG_FATAL(message_expr) ALGLIB_LOG_STREAM(::alglib::logging::Level::Fatal, message_expr)

