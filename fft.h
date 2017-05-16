#ifndef _FFT_H
#define _FFT_H

#include "macros.h"

typedef struct _COMPLEX {
	double real;
	double imag;
} COMPLEX, *PCOMPLEX;



COMPLEX* _cdecl fft(COMPLEX* f, uint32_t n, int32_t s);
COMPLEX* realfft(COMPLEX* f, uint32_t n, int32_t s);

#endif