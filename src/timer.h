/**
 * @file
 *
 * Simple timer class
 */

#ifndef VR_TIMER_H_
#define VR_TIMER_H_

#include <ostream>
#include <chrono>

#include "ioutils.h"


namespace vr {

/**
 * A simple timer class that starts measuring time when created and returns
 * the elapsed time when requested.
 */
class Timer {

public:

	using clock = std::chrono::high_resolution_clock;
	using duration = typename std::chrono::microseconds::rep;

	/**
	 * Returns the number of milliseconds elapsed since the creation
	 * of the timer
	 *
	 * @return The time elapsed since the creation of the timer, in [us]
	 */
	inline
	duration get() const {
		return std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - t0).count();
	}

private:
	clock::time_point t0 {clock::now()};

};

template <typename T>
inline
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const Timer &t) {
	os << us_time(t.get());
	return os;
}

}  // namespace vr

#endif // VR_TIMER_H_