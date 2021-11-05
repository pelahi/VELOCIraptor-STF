#ifndef VR_LOGGING_H_
#define VR_LOGGING_H_

#include <iosfwd>
#include <string>

#include "allvars.h"

namespace vr
{

/// Supported logging levels
enum LogLevel : char {
	trace = 0,
	debug,
	info,
	warning,
	error
};

/// Convert the given LogLevel into a string
std::string to_string(LogLevel level);

/// Initializes the logging system
void init_logging(LogLevel log_level);

/// Checks if the given log level would result on a statement being printed
bool log_enabled(LogLevel log_level);

/// A class that streams log statements to stdout
class LogStatement {

private:
	std::ostream &m_os;
	LogLevel m_level;
	bool m_enabled = true;

public:
	LogStatement(LogLevel level, const char *file, int line);
	LogStatement();
	~LogStatement() noexcept;

	template <typename T>
	LogStatement &operator<<(const T &v)
	{
		if (m_enabled && log_enabled(m_level)) {
			m_os << v;
		}
		return *this;
	}
};

}  // namespace vr

/// Macro to perform logging via output stream operator
#define LOG(lvl) vr::LogStatement(vr::LogLevel::lvl, __FILE__, __LINE__)

/// Macro to perform logging via output stream operator on the first MPI rank only
#ifdef USEMPI
#define LOG_RANK0(lvl) ((ThisTask == 0) ? LOG(lvl) : vr::LogStatement())
#else
#define LOG_RANK0(lvl) LOG(lvl)
#endif // USEMPI

/// Macro used to check if a logging level is enabled, used only when logging
/// is really expensive
#define LOG_ENABLED(lvl) log_enabled(vr::LogLevel::lvl)

#endif /* VR_LOGGING_H_ */