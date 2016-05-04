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

#ifdef __cplusplus
#ifndef typeof
#define typeof(x) decltype(x)
#endif
#endif
