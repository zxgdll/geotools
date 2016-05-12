#ifndef __GEOTOOLS_H__
#define __GEOTOOLS_H__

#include <limits>

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

#define _sec2rad(x) (x * 4.84813681 / 1000000000.0)
#define _deg(x) (x * 180.0 / PI)
#define _rad(x) (x * PI / 180.0)

#endif
