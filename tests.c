#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "bint.h"


word randw() {
    word w = rand();
    while (!(w >> (WORDSIZE - 2))) w *= rand();
    return w;
}

void randomize_bigint(bint *a) {
    a->len = rand() % a->mem;
    for (ulong k = 0; k < a->len; k++) {
        while (!a->digits[k])
            a->digits[k] = randw();
    }
    a->sgn = ((ulong)rand()) % 2;
}

void fail(bint *s, bint *t, bint *a, bint *b, bint *c, const char *error) {
    printf("ERROR: %s:\n", error);
    printf("L = %s\n", bint_to_str(s, 10));
    printf("R = %s\n", bint_to_str(t, 10));
    printf("SEED:\n");
    printf("a = %s\n", bint_to_str(a, 10));
    printf("b = %s\n", bint_to_str(b, 10));
    printf("c = %s\n", bint_to_str(c, 10));
}


bool test_additive_laws(bmc *ctx, ulong rounds) {
    bool sane = true;
    sint shorty = 0;

    bint *a = bint_init_alloc(ctx, 0, 9);
    bint *b = bint_init_alloc(ctx, 0, 9);
    bint *c = bint_init_alloc(ctx, 0, 9);
    bint *t = bint_init(ctx, 0),
        *r = bint_init(ctx, 0),
        *s = bint_init(ctx, 0);
    bint *vars[6] = { s,r,t,c,b,a };

    for (ulong k = 0; k < rounds; k++) {
        randomize_bigint(a);
        randomize_bigint(b);
        randomize_bigint(c);

        bint_add(t, a, b);
        bint_sub(s, t, b);
        if (bint_compare(a, s) != 0) {
            fail(s, a, a, b, c, "(a+b)-b == a");
            sane = false;
            break;
        }

        bint_mul_to(t, c);
        bint_mul(r, a, c);
        bint_mul(s, c, b);
        bint_add_to(s, r);
        if (bint_compare(t, s) != 0) {
            fail(t, s, a, b, c, "(a+b)*c == (a*c)+(c*b)");
            sane = false;
            break;
        }

        if (!bint_is_zero(b)) {
            b->sgn = 0;
            bint_div(t, c, b);
            bint_mod(s, c, b);
            bint_mul(r, t, b);
            bint_add_to(r, s);
            if (bint_compare(r, c) != 0) {
                fail(r, c, a, b, c, "b*(c/b) + (c%b) == c");
                sane = false;
                break;
            }
        }

        while (!shorty)
            shorty = (sint)randw();
        if (shorty < 0) shorty = -shorty;

        bint_short_assign(s, shorty);
        bint_mod(r, a, s);
        bint_mod(t, b, s);
        bint_add_to(r, t);
        bint_mod_to(r, s);
        bint_add(t, a, b);
        bint_mod_to(t, s);
        if (bint_compare(t, r) != 0) {
            fail(r, t, a, b, s, "((a%c)+(b%c))%c == (a+b)%c");
            sane = false;
            break;
        }
    }


    for (int i = 0; i < 6; i++)
        free(vars[i]);

    return sane;
}

bool report(bool passed, const char *what) {
    printf("[%s] %s\n", passed ? "pass" : "fail", what);
    return passed;
}

bool test(bmc *ctx) {
    bool passed = true;
    passed = report(test_additive_laws(ctx, 20000), "standard laws of arithmetic");

    putchar('\n');
    return passed;
}

int main(int argc, char** argv) {
    srand((unsigned int)time(NULL));
    bmc *ctx = bmc_create();

    printf("[info] digits have %d bits.\n", WORDSIZE);

    printf("[info] performing tests for standard settings.\n");
    test(ctx);

    printf("[info] performing tests for FFT-only settings.\n");
    test(ctx);

    getchar();
    bmc_destroy(ctx);
    return 0;
}