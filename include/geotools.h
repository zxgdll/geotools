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

#ifdef _MSC_VER
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

#define G_PI 3.14159265358979323846
#define G_E 2.71828

#define G_DBL_MAX_POS (std::numeric_limits<double>::max())
#define G_DBL_MAX_NEG (std::numeric_limits<double>::lowest())

#define G_FLT_MAX_POS (std::numeric_limits<float>::max())
#define G_FLT_MAX_NEG (std::numeric_limits<float>::lowest())

#define g_min(a, b) ((a) > (b) ? (b) : (a))
#define g_max(a, b) ((a) < (b) ? (b) : (a))
#define g_sq(a) ((a) * (a))
#define g_abs(x) ((x) < 0 ? -(x) : (x))
#define g_deg(x) ((x) * 180.0 / G_PI)
#define g_rad(x) ((x) * G_PI / 180.0)

extern int g__loglevel;

#define G_LOG_TRACE 5
#define G_LOG_DEBUG 4
#define G_LOG_WARN 3
#define G_LOG_ERROR 2
#define G_LOG_NONE 0

#define g_loglevel(x) {g__loglevel = x;}
/*
#define g_log(x, y) _Pragma("omp critical") \
        { \
                if(g__loglevel >= y) { \
                        std::cerr << std::setprecision(12) << x << std::endl; \
                } \
        }
 */
#define g_log(x, y) { if(g__loglevel >= y) std::cerr << std::setprecision(12) << x << std::endl; }
#define g_trace(x) g_log("TRACE:   " << x, G_LOG_TRACE)
#define g_debug(x) g_log("DEBUG:   " << x, G_LOG_DEBUG)
#define g_warn(x)  g_log("WARNING: " << x, G_LOG_WARN)
#define g_error(x) g_log("ERROR:   " << x, G_LOG_ERROR)

#define g_argerr(x) {std::stringstream _ss; _ss << x; throw std::invalid_argument(_ss.str());}
#define g_implerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}
#define g_runerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}

#endif
