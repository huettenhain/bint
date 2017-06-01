#include "fft.h"
#include <math.h>

COMPLEX* realfft(COMPLEX* g, uint32_t n, int32_t s) {

    unsigned long k, N = n >> 1;
    double tmpval;
    double a, b, c, d;
    COMPLEX w, e;

    s /= abs(s);

    tmpval = s * 3.1415926535897931 / n;
    e.real = w.real = cos(tmpval);
    e.imag = w.imag = sin(tmpval);

    if (s < 0) {
        fft(g, n, s);

        g[0].real /= 2;
        g[0].imag /= 2;
        g[N].real /= 2;
        g[N].imag = -g[N].imag / 2;

        for (k = 1; k < N; k++) {
            g[k].real /= 4;
            g[k].imag /= 4;
        }
        for (k = n - 1; k > N; k--) {
            g[k].real /= 4;
            g[k].imag /= 4;
        }
    }
    else {
        g[N].real = 2 * g[N].real;
        g[N].imag = -2 * g[N].imag;
    }

    for (s = -s, k = 1; k < N; k++)
    {
        a = (g[k].real + g[n - k].real);
        b = (g[k].real - g[n - k].real) * s;
        c = (g[k].imag + g[n - k].imag) * s;
        d = (g[k].imag - g[n - k].imag);

        g[k].real = a + b*e.imag + c*e.real;
        g[k].imag = d - b*e.real + c*e.imag;
        g[n - k].real = a - b*e.imag - c*e.real;
        g[n - k].imag = c*e.imag - b*e.real - d;

        e.real = (tmpval = e.real)*w.real - e.imag*w.imag;
        e.imag = e.imag *w.real + tmpval*w.imag;
    }

    g[0].real = (tmpval = g[0].real) + g[0].imag;
    g[0].imag = tmpval - g[0].imag;

    return (s < 0) ? fft(g, n, 1) : g;
}

