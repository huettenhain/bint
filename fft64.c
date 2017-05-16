#include "fft.h"

#ifdef BITS64

#include <stdlib.h>
#include <memory.h>
#include <malloc.h>
#include <math.h>

COMPLEX* fft(COMPLEX* f, uint32_t n, int32_t s) {

	uint32_t k, j, M, N = n;
	double tmpval;
	COMPLEX w, z, a;

	s /= abs(s);
	tmpval = s * 6.2831853071795862 / n;

	w.real = cos(tmpval);
	w.imag = sin(tmpval);

	if (s < 0) for (j = 0; j<n; j++) {
		f[j].imag /= n;
		f[j].real /= n;
	}

	while (N > 1) {

		M = N >> 1;

		for (j = 0; j<n; j += N)
		{
			z.real = 1;
			z.imag = 0;

			for (k = 0; k<M; k++)
			{
				a = f[j + k];
				a.real -= f[j + k + M].real;
				a.imag -= f[j + k + M].imag;
				tmpval = a.real;
				a.real = tmpval*z.real - a.imag*z.imag;
				a.imag = a.imag*z.real + tmpval*z.imag;

				f[j + k].real += f[j + k + M].real;
				f[j + k].imag += f[j + k + M].imag;

				f[j + k + M] = a;

				tmpval = z.real;
				z.real = z.real*w.real - z.imag*w.imag;
				z.imag = z.imag*w.real + tmpval*w.imag;
			}
		}

		tmpval = w.real;
		w.real = tmpval*w.real - w.imag*w.imag;
		w.imag *= tmpval;
		w.imag += w.imag;

		N = M;
	}


	for (j = k = 0; k<n; k++) {
		if (j<k) {
			a = f[j];
			f[j] = f[k];
			f[k] = a;
		}
		for (N = n; j<N; j ^= N)
			N >>= 1;
	}

	return f;
}

#endif 
