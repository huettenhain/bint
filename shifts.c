#include <string.h>
#include "bint.h"

ulong braw_lshift(word *a, ulong size, ulong places) {
    ulong shiftOffset = places % WORDSIZE;
    ulong otherShift = WORDSIZE - shiftOffset;
    ulong delta = places / WORDSIZE;
    ulong words_null = delta;

    if (delta < size) {
        if (shiftOffset) {
            word thisWord = a[size - delta - 1];
            ulong i, wp = size - 1;
            for (i = size - delta - 2; i < size; i--) {
                a[wp--] = (thisWord << shiftOffset) | (a[i] >> otherShift);
                thisWord = a[i];
            }
            if ((a[wp] = thisWord << shiftOffset) == 0 && wp)
                words_null++;
        }
        else {
            memmove(a + delta, a, (size - delta) * sizeof(word));
        }
        memset(a, 0, delta * sizeof(word));
        return words_null;
    }
    else {
        memset(a, 0, size * sizeof(word));
        return size - 1;
    }
}

ulong braw_rshift(word *a, ulong size, ulong places) {
    ulong shiftOffset = places % WORDSIZE;
    ulong otherShift = WORDSIZE - shiftOffset;
    ulong delta = places / WORDSIZE;
    ulong words_null = delta;

    if (delta < size) {
        if (shiftOffset) {
            word thisWord = a[delta];
            ulong i, wp = 0;
            for (i = delta + 1; i < size; i++) {
                a[wp++] = (thisWord >> shiftOffset) | (a[i] << otherShift);
                thisWord = a[i];
            }
            if ((a[wp] = thisWord >> shiftOffset) == 0 && wp)
                words_null++;
        }
        else {
            memmove(a, a + delta, (size - delta) * sizeof(word));
        }
        memset(a + size - delta, 0, delta * sizeof(word));
        return words_null;
    }
    else {
        memset(a, 0, size * sizeof(word));
        return size - 1;
    }
}

void braw_rshift1(word *a, ulong size) {
    word *p;
    ASSERT(a);
    p = a + size - 1;
    while (a < p) {
        *a >>= 1; a++;
        *(a - 1) |= *a << (WORDSIZE - 1);
    }
    *a >>= 1;
}

void braw_lshift1(word *a, ulong size) {
    word *p;
    ASSERT(a);
    p = a + size - 1;
    while (p > a) {
        *p <<= 1; p--;
        *(p + 1) |= *p >> (WORDSIZE - 1);
    }
    *p <<= 1;
}

bint *bint_lshift(bint *a, ulong places) {
    ulong needed_bits;
    ASSERT(a);
    needed_bits = places - WORDSIZE + word_log2(a->digits[a->len - 1]) + 1;
    if (needed_bits <= places) {
        ulong needed_words = needed_bits / WORDSIZE;
        if (needed_bits % WORDSIZE) needed_words++;
        if (!bint_alloc(a, a->len + needed_words))
            return NULL;
        else a->len += needed_words;
    }
    braw_lshift(a->digits, a->len, places);
    return a;
}

bint *bint_rshift(bint *a, ulong places) {
    ASSERT(a);
    a->len -= braw_rshift(a->digits, a->len, places);
    return a;
}

bint *bint_lshift1(bint *a) {
    if ((a->digits[a->len - 1] & BINT_MSB_MASK) == BINT_MSB_MASK) {
        if (!bint_alloc(a, a->len + 1)) return NULL;
        else a->len++;
    }

    braw_lshift1(a->digits, a->len);
    return a;
}

bint *bint_rshift1(bint *a) {
    braw_rshift1(a->digits, a->len);
    if (!a->digits[a->len - 1]) a->len--;
    return a;
}
