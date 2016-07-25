#ifndef __GEOTOOLS_H__
#define __GEOTOOLS_H__

#include <limits>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>

#ifdef __GNUC__
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

#define PI 3.14159265358979323846
#define E 2.71828

#define DBL_MAX_POS (std::numeric_limits<double>::max())
#define DBL_MAX_NEG (std::numeric_limits<double>::lowest())

#define FLT_MAX_POS (std::numeric_limits<float>::max())
#define FLT_MAX_NEG (std::numeric_limits<float>::lowest())

#define _min(a, b) ((a) > (b) ? (b) : (a))
#define _max(a, b) ((a) < (b) ? (b) : (a))
#define _sq(a) ((a) * (a))
#define _abs(x) (x < 0 ? -x : x)
#define _deg(x) (x * 180.0 / PI)
#define _rad(x) (x * PI / 180.0)

static int __loglevel = 0;

#define LOG_TRACE 5
#define LOG_DEBUG 4
#define LOG_WARN 3
#define LOG_ERROR 2
#define LOG_NONE 0

#define _loglevel(x) {__loglevel = x;}
#define _log(x, y) if(__loglevel >= y) { std::cerr << x << std::endl; }
#define _trace(x) _log("TRACE:   " << x, LOG_TRACE)
#define _debug(x) _log("DEBUG:   " << x, LOG_DEBUG)
#define _warn(x)  _log("WARNING: " << x, LOG_WARN)
#define _error(x) _log("ERROR:   " << x, LOG_ERROR)

#define _argerr(x) {std::stringstream _ss; _ss << x; throw std::invalid_argument(_ss.str());}
#define _implerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}
#define _runerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}

static std::string __status_buf;

void _status(int step, int of) {
	if(step > of)
		of = step;
	for(int i = 0; i < __status_buf.size(); ++i)
		std::cerr << '\b';
	std::stringstream buf;
	buf << "Status: " << std::setprecision(2) << (of <= 0 ? 0 : (float) step / of) << "%%" << std::endl;
	__status_buf.assign(buf.str());
	std::cerr << __status_buf;
}
#endif
