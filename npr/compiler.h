#ifndef NPR_COMPILER_HPP
#define NPR_COMPILER_HPP

#ifdef __GNUC__
#define ALWAYS_INLINE __attribute__((always_inline)) __inline
#elif defined _MSC_VER
#define ALWAYS_INLINE __forceinline
#else
#define ALWAYS_INLINE __inline
#endif
#endif
