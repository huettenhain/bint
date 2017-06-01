#include <malloc.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "bint.h"

void randomize_buffer(word *b, ulong l) {
    ulong k;
    for (; l; l--) {
        b[l - 1] = rand();
        for (k = 0; k < WORDSIZE; k += 16)
            b[l - 1] *= rand();
    }
}

bool bmc_resize(bmc *context, ulong size) {
    void *tmp;
    if (context->size >= size)
        return true;
    if (tmp = realloc(context->buffer, size)) {
        context->buffer = tmp;
        context->size = size;
        return true;
    }
    else {
        return false;
    }
}

bmc *bmc_create() {
    bmc *context = malloc(sizeof(bmc));
    if (context) {
        ulong len = KMUL_BUFFERSIZE(BINT_INITSIZE);
        if (context->buffer = malloc(len)) {
            context->size = len;
        }
        else {
            free(context);
            context = NULL;
        }
    }
    bmc_defaults(context);
    return context;
}

void bmc_defaults(bmc *context) {
    ASSERT(context);
    ulong i;
    context->k_cut = BINT_KCUT;
#   define FFT_IGNOR 6
#   define FFT_TIMID 10
#   define FFT_EAGER 15
#   define FFT_BLIND 20
    for (i = 0; i < FFT_IGNOR; i++)
        context->f_dst[i] = -1;
    for (i = FFT_IGNOR; i < FFT_TIMID; i++) context->f_dst[i] = TWO_TO_THE(i - 4);
    for (i = FFT_TIMID; i < FFT_EAGER; i++) context->f_dst[i] = TWO_TO_THE(i - 3);
    for (i = FFT_EAGER; i < FFT_BLIND; i++) context->f_dst[i] = TWO_TO_THE(i - 2);
    for (i = FFT_BLIND; i < WORDSIZE - 1; i++) context->f_dst[i] = TWO_TO_THE(i - 1);
}
void bmc_fft(bmc *context) {
    ulong i;
    for (i = 0; i < WORDSIZE; i++)
        context->f_dst[i] = TWO_TO_THE(WORDSIZE - 2);
}

void bmc_destroy(bmc *context) {
    if (context->buffer)
        free(context->buffer);
    free(context);
}