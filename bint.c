#include <malloc.h>
#include <memory.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "bint.h"

bint *bint_init(bmc *context, sint init) {
    return bint_init_alloc(context, init, BINT_INITSIZE);
}

void  bint_fit(bint *a) {
    while (a->len && !a->digits[a->len - 1])
        a->len--;
}

bint *bint_init_alloc(bmc *context, sint init, ulong size) {
    bint *a;

    if (!context || !context->buffer || !(a = malloc(sizeof(bint))))
        return NULL;

    a->len = 0;
    a->mem = 0;
    a->sgn = 0;
    a->context = context;

    if (bint_alloc(a, size)) {
        if (init < 0) {
            a->sgn = 1;
            init = -init;
        }
        a->digits[0] = (word)init;
        a->len = 1;
        return a;
    }
    else {
        free(a);
        return NULL;
    }
}

char *bint_to_str(bint *a, ulong base) {
    word digit, pos = 0;
    bint *t;
    char *s = NULL;

    if (bint_is_zero(a)) {
        return "0";
    }

    if ((t = bint_init(a->context, 0))) {
        if (bint_assign(t, a)) {
            ulong bytes = (ulong)ceil((bint_log2(a) + 1) * log(2) / log(base)) + 2;
            if (bmc_resize(a->context, bytes)) {
                s = (char*)a->context->buffer;
                while (!bint_is_zero(t)) {
                    bint_short_divu_to(t, base, &digit);
                    s[pos++] = (char)(digit < 10 ?
                        ('0' + digit) : ('A' + digit - 10));
                }
                if (t->sgn) s[pos++] = '-';
                s[pos] = '\0';
                for (word k = 0; k < pos / 2; k++) {
                    char tmp = s[k];
                    s[k] = s[pos - k - 1];
                    s[pos - k - 1] = tmp;
                }
            }
        }
        bint_free(t);
    }

    return s;
}

bint *bint_from_str(bint *a, char *s) {

    ulong l = (ulong)strlen(s);
    ulong sgn = 0;
    if (!l) return NULL;

    if (*s == '-')
        sgn = 1, s++, l--;

    if (l < 2) {
        if (!l || !isdigit(*s))
            return NULL;
        else return bint_short_assign(a, *s - '0');
    }
    else {
        word  base = 10;
        ulong size = (ulong)ceil(3.3219280948873626*l);
        char  b[2] = { '\0', '\0' };
        if (*s == '0') {
            switch (s[1]) {
            case'b': base = 0x02;
                size = l;
                s += 2;
                break;
            case'o': base = 0x08;
                size = 3 * l;
                s += 2;
                break;
            case'x': base = 0x10;
                size = 4 * l;
                s += 2;
                break;
            }
        }

        size = size % WORDSIZE ? (size / WORDSIZE) + 1 :
            (size / WORDSIZE);

        bint_short_assign(a, 0);

        if (!bint_alloc(a, size))
            return NULL;

        else for (; *s; s++) {
            b[0] = *s;
            bint_short_mulu_to(a, base);
            bint_short_add_to(a, strtol(b, NULL, (int)base));
        }

        a->sgn = sgn;
        return a;
    }
}

bint *bint_alloc(bint *a, ulong size) {

    ASSERT(a);

    if (a->mem >= size)
        return a;

    if (!a->mem) {
        if (a->digits = calloc(size, sizeof(word))) {
            a->mem = size;
            return a;
        }
        else {
            a->mem = 0;
            return NULL;
        }
    }
    else {
        word *tmp;
        if (tmp = realloc(a->digits, sizeof(word)*size)) {
            memset(
                (a->digits = tmp) + a->mem, 0,
                (size - a->mem) * sizeof(word));
            a->mem = size;
            return a;
        }
        else {
            return NULL;
        }
    }
}


void bint_free(bint *a) {
    ASSERT(a);
    if (a->digits) {
        free(a->digits);
        a->digits = NULL;
    }
    free(a);
}

ulong word_log2(word a) {
    unsigned long l = 0 - 1;
    do l++; while (a >>= 1);
    return l;
}

ulong braw_logB(word *x, ulong n) {
    for (--n; n; n--) if (x[n]) return n;
    return 0;
}

ulong braw_log2(word *x, ulong n) {
    ulong k = braw_logB(x, n);
    return (k*WORDSIZE) + word_log2(x[k]);
}

ulong bint_logB(bint *a) {
    return braw_logB(a->digits, a->len);
}

ulong bint_log2(bint *a) {
    return braw_log2(a->digits, a->len);
}

ulong word_log2_ceil(word a) {
    ulong p = word_log2(a);
    if (a > (((word)1) << p)) p++;
    return p;
}

word word_next_pow2(word a) {
    ulong p = 1 << word_log2(a);
    if (a > p) p <<= 1;
    return p;
}

bool bint_is_zero(bint *a) {
    ASSERT(a);
    return (a->len == 0);
}

bool bint_is_sint(bint *a) {
    ASSERT(a);
    return (a->len == 1) && (a->digits[0] < TWO_TO_THE(WORDSIZE - 1));
}
sint bint_to_sint(bint *a) {
    ASSERT(a);
    ASSERT(bint_is_sint(a));
    sint r = (sint)a->digits[0];
    if (a->sgn) return -r;
    else return r;
}

int braw_compare(word *a, word *b, ulong n, ulong m) {
    if (n < m) return -1;
    else if (n > m) return  1;
    while (n--) {
        if (a[n] < b[n]) return -1;
        else if (a[n] > b[n]) return  1;
    }
    return 0;
}

int bint_compare(bint *a, bint *b) {
    if (bint_is_zero(a)) {
        if (bint_is_zero(b)) return 0;
        else if (b->sgn) return 1;
        else return -1;
    }
    else if (bint_is_zero(b)) {
        if (a->sgn) return -1;
        else return 1;
    }

    if (a->sgn ^ b->sgn)
        return (2 * ((sint)b->sgn)) - 1;
    else return b->sgn ?
        braw_compare(b->digits, a->digits, b->len, a->len)
        : braw_compare(a->digits, b->digits, a->len, b->len);
}
