#ifndef __GEOTOOLS_H__
#define __GEOTOOLS_H__

#ifdef __GNUC__
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

#define PI 3.14159265358979323846
#define _min(a, b) (a > b ? b : a)
#define _max(a, b) (a < b ? b : a)
#define _sq(a) (a * a)
#endif
