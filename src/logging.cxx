#include <iomanip>
#include <stdexcept>
#include <iostream>

#include "allvars.h"
#include "proto.h"
#include "logging.h"
#include "timer.h"

namespace vr
{

/// The log level with which the library is configured
static LogLevel configured_level {LogLevel::info};

/// A timer that starts at the beginning of the program execution
static Timer start_timer;

std::string to_string(LogLevel level)
{
	switch (level) {
	case trace:
		return "trace";
	case debug:
		return "debug";
	case info:
		return "info";
	case warning:
		return "warn";
	case error:
		return "error";
	default:
		throw std::invalid_argument("Unsupported level: " + std::to_string(static_cast<char>(level)));
	}
}

bool log_enabled(LogLevel level)
{
	return level >= configured_level;
}

void init_logging(LogLevel lvl)
{
	configured_level = lvl;
}

LogStatement::LogStatement()
  : m_os(std::cout), m_level(LogLevel::error), m_enabled(false)
{
}

LogStatement::LogStatement(LogLevel level, const char *file, int line)
  : m_os(std::cout), m_level(level)
{
	if (!log_enabled(m_level)) {
		return;
	}
#ifdef USEMPI
	m_os << '[' << std::setw(4) << std::setfill('0') << ThisTask << "] ";
#endif
	auto us_since_start = start_timer.get();
	auto secs = us_since_start / 1000000;
	auto msecs = (us_since_start - (secs * 1000000)) / 1000;
	m_os << '['
	   << std::setw(4) << std::setfill(' ') << secs << '.'
	   << std::setw(3) << std::setfill('0') << msecs << "] ["
	   << std::setw(5) << std::setfill(' ') << to_string(m_level) << "] ";
#ifdef VR_LOG_SOURCE_LOCATION
	m_os << basename(file) << ':' << line << ' ';
#endif // VR_LOG_SOURCE_LOCATION
}

LogStatement::~LogStatement() noexcept
{
	if (!m_enabled || !log_enabled(m_level)) {
		return;
	}
	try {
		m_os << '\n';
	} catch (const std::exception &) {
		// silently ignore...
	}
}

}  // namespace vr