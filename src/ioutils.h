/**
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia
 * Copyright by UWA(in the framework of the ICRAR)
 * All rights reserved
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111 - 1307  USA
 */

/**
 * @file
 *
 * Utilities for uniformly formatting values
 */

#ifndef VR_IOUTILS_H
#define VR_IOUTILS_H

#include <cassert>
#include <chrono>
#include <iomanip>
#include <ostream>
#include <vector>

namespace vr
{

namespace detail {

	template <int N, typename T>
	struct _fixed {
		T _val;
	};

	template <typename T, int N, typename VT>
	inline
	std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, detail::_fixed<N, VT> v)
	{
		os << std::setprecision(N) << std::fixed << v._val;
		return os;
	}

} // namespace detail

///
/// Sent to a stream object, this manipulator will print the given value with a
/// precision of N decimal places.
///
/// @param v The value to send to the stream
///
template <int N, typename T>
inline
detail::_fixed<N, T> fixed(T v) {
	return {v};
}

namespace detail {

	struct _memory_amount {
		std::size_t _val;
	};

	struct _microseconds_amount {
		std::chrono::microseconds::rep _val;
	};

	template <typename ForwardIterator>
	struct _printable_range {
		ForwardIterator begin;
		ForwardIterator end;
		std::string sep;
	};

	template <typename T>
	inline
	std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const detail::_memory_amount &m)
	{

		if (m._val < 1024) {
			os << m._val << " [B]";
			return os;
		}

		float v = m._val / 1024.;
		const char *suffix = " [KiB]";

		if (v > 1024) {
			v /= 1024;
			suffix = " [MiB]";
		}
		if (v > 1024) {
			v /= 1024;
			suffix = " [GiB]";
		}
		if (v > 1024) {
			v /= 1024;
			suffix = " [TiB]";
		}
		if (v > 1024) {
			v /= 1024;
			suffix = " [PiB]";
		}
		if (v > 1024) {
			v /= 1024;
			suffix = " [EiB]";
		}
		// that should be enough...

		os << fixed<3>(v) << suffix;
		return os;
	}

	template <typename T>
	inline
	std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const detail::_microseconds_amount &t)
	{
		auto time = t._val;
		if (time < 1000) {
			os << time << " [us]";
			return os;
		}

		time /= 1000;
		if (time < 1000) {
			os << time << " [ms]";
			return os;
		}

		float ftime = time / 1000.f;
		const char *prefix = " [s]";
		if (ftime > 60) {
			ftime /= 60;
			prefix = " [min]";
			if (ftime > 60) {
				ftime /= 60;
				prefix = " [h]";
				if (ftime > 24) {
					ftime /= 24;
					prefix = " [d]";
				}
			}
		}
		// that should be enough...

		os << fixed<3>(ftime) << prefix;
		return os;
	}

	template <typename T, typename ForwardIterator>
	std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os,
		const _printable_range<ForwardIterator> &range)
	{
		auto first = range.begin;
		os << '[';
		for (; first != range.end;) {
			os << *first;
			first++;
			if (first != range.end) {
				os << range.sep;
			}
		}
		os << ']';
		return os;
	}

} // namespace detail

///
/// Sent to a stream object, this manipulator will print the given amount of
/// memory using the correct suffix and 3 decimal places.
///
/// @param v The value to send to the stream
///
inline
detail::_memory_amount memory_amount(std::size_t amount) {
	return {amount};
}

///
/// Sent to a stream object, this manipulator will print the given amount of
/// nanoseconds using the correct suffix and 3 decimal places.
///
/// @param v The value to send to the stream
///
inline
detail::_microseconds_amount us_time(std::chrono::microseconds::rep amount) {
	return {amount};
}

///
/// Sent to a stream object, this manipulator will print the given range of
/// values into
///
/// @param v The value to send to the stream
///
template <typename Value>
detail::_printable_range<typename std::vector<Value>::const_iterator>
printable_range(const std::vector<Value> &values, const std::string &sep = ", ")
{
	return {values.cbegin(), values.cend(), sep};
}

template <typename Value>
detail::_printable_range<const Value *>
printable_range(const Value *values, std::size_t size, const std::string &sep = ", ")
{
	return {values, values + size, sep};
}

}  // namespace vr

#endif // VR_IOUTILS_H
