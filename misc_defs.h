#ifdef _WIN32
#ifdef _WIN64
typedef unsigned __int64 size_t;
typedef signed __int64 ssize_t;
#else
typedef _W64 unsigned int size_t;
typedef _W64 signed int ssize_t;
#endif
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962
#endif

#ifdef __cplusplus
#ifndef typeof
#define typeof(x) decltype(x)
#endif
#endif
