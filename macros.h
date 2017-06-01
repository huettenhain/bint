#if                                                    0
#
# This is where we try to figure out which architecture
# and compiler we are on, so we can include optimized 
# ASM instructions later in the code. Since even on the
# same architecture, inline ASM has to differ from
# compiler to compiler, it is a big TODO to add more of
# these.
#
# Also, some useful macros are defined.
#                                                  endif
#
#ifndef _MACROS_H
#define _MACROS_H
#
#undef BITS32
#undef BITS64
#
#if defined(_MSC_VER)
# if defined(_M_IX86) 
#  define MVC_X86_32
#  define BITS32
# elif defined(_M_X64)
#  define MVC_X86_64
#  define BITS64
# endif
#elif defined(__GNUC__)
# if defined(__i386__)
#  define GCC_X86_32
#  define BITS32
# elif defined(__x86_64__)
#  define GCC_X86_64
#  define BITS64
# endif
#endif
#
# ifndef DEBUG
#  define _NDEBUG
# endif
#
# include <assert.h>
#
# ifdef DEBUG
#  ifdef _WIN32
#   define ASSERT(_expr) __assume(_expr)
#  else
#   define ASSERT(_expr) (void)0
#  endif
# else
#  define ASSERT(_expr) assert(_expr)
# endif
#
# ifdef DEBUG
#  define NODEFAULT default: ASSERT(0)
# else
#  ifdef _WIN32
#   define NODEFAULT default: __assume(0)
#  else
#   define NODEFAULT
#  endif
# endif
#
# ifdef _WIN32
#  define INLINE __inline
# else
#  define INLINE inline
# endif
#
# define MAX(_x,_y)   ( ((_x)>(_y)) ? (_x) : (_y) )
# define MIN(_x,_y)   ( ((_x)<(_y)) ? (_x) : (_y) )
# define DIVUP(_x,_y) ( (_x)%(_y) ? (_x)/(_y)+1 : (_x)/(_y) )
# 
# if defined(BITS64)
#   define TWO_TO_THE(_x) ( 1i64 << (_x) )
# elif defined(BITS32)
#   define TWO_TO_THE(_x) ( 1i32 << (_x) )
# else 
#   define TWO_TO_THE(_x) ( 1    << (_x) )
# endif 
#
#ifdef _MSC_VER
 typedef   signed __int64  int64_t;
 typedef   signed __int32  int32_t;
 typedef   signed __int16  int16_t;
 typedef   signed __int8    int8_t;
 typedef unsigned __int64 uint64_t;
 typedef unsigned __int32 uint32_t;
 typedef unsigned __int16 uint16_t;
 typedef unsigned __int8   uint8_t;
#else
# include <stdint.h>
#endif
#endif
