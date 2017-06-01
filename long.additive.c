#include <string.h>

#include "bint.h"
#include "atomics_inline.c"

word braw_shorter_add(word *dst, word *a, word *b, ulong a_len, ulong b_len) {
    word carry;

    ASSERT(dst && a && b);
    ASSERT(b_len <= a_len);

    a_len -= b_len;

    for (carry = 0; b_len; b_len--)
        *dst++ = add_with_carry(*a++, *b++, &carry);
    for (; carry && a_len; a_len--)
        *dst++ = add_with_carry(*a++, 0, &carry);

    memcpy(dst, a, a_len * sizeof(word));
    return carry;
}

word braw_shorter_add_to(word *a, word *b, ulong a_len, ulong b_len) {
    word carry;

    ASSERT(a && b);
    ASSERT(b_len <= a_len);

    a_len -= b_len;
    for (carry = 0; b_len; b_len--, a++)
        *a = add_with_carry(*a, *b++, &carry);
    for (; carry && a_len; a_len--, a++)
        *a = add_with_carry(*a, 0, &carry);

    return carry;
}

word braw_shorter_sub(word *dst, word *a, word *b, ulong a_len, ulong b_len) {
    word borrow;

    ASSERT(dst && a && b);
    ASSERT(b_len <= a_len);

    a_len -= b_len;
    for (borrow = 0; b_len; b_len--)
        *dst++ = sub_with_borrow(*a++, *b++, &borrow);
    for (; borrow && a_len; a_len--)
        *dst++ = sub_with_borrow(*a++, 0, &borrow);

    memcpy(dst, a, a_len * sizeof(word));
    return borrow;
}

word braw_shorter_sub_to(word *a, word *b, ulong a_len, ulong b_len) {
    word borrow;

    ASSERT(a && b);
    ASSERT(b_len <= a_len);

    a_len -= b_len;
    for (borrow = 0; b_len; b_len--, a++)
        *a = sub_with_borrow(*a, *b++, &borrow);
    for (; borrow && a_len; a_len--, a++)
        *a = sub_with_borrow(*a, 0, &borrow);

    return borrow;
}

word braw_add(word *dst, word *a, word *b, ulong len) {
    word carry;
    ASSERT(dst && a && b);
    for (carry = 0; len; len--)
        *dst++ = add_with_carry(*a++, *b++, &carry);
    return carry;
}

word braw_add_to(word *a, word *b, ulong len) {
    return braw_add(a, a, b, len);
}

word braw_sub(word *dst, word *a, word *b, ulong len) {
    word borrow;
    ASSERT(dst && a && b);
    for (borrow = 0; len; len--)
        *dst++ = sub_with_borrow(*a++, *b++, &borrow);
    return borrow;
}

word braw_sub_to(word *a, word *b, ulong len) {
    return braw_sub(a, a, b, len);
}

bint *_bint_add_digits(bint *result, bint *a, bint *b, ulong sub) {
    ASSERT(result && a && b);
    word carry;
    ulong sgn = 0;

    if (a->len < b->len)
        return _bint_add_digits(result, b, a, sub);

    if (sub) {
        if (!bint_alloc(result, a->len))
            return NULL;
        else {
            carry = braw_shorter_sub(result->digits,
                a->digits, b->digits, a->len, b->len);
            result->len = a->len;
            if (carry) {
                braw_neg(result->digits, a->len);
                sgn = 1;
            }
        }
    }
    else {
        if (!bint_alloc(result, a->len + 1))
            return NULL;
        else {
            carry = braw_shorter_add(result->digits,
                a->digits, b->digits, a->len, b->len);
            result->len = a->len;
            if (carry) {
                result->digits[a->len] = 1; // carry has to be 1
                result->len++;
            }
        }
    }
    bint_fit(result);
    result->sgn = a->sgn ^ sgn;
    return result;
}

bint *bint_assign(bint *a, bint *b) {
    ASSERT(a && b);

    if (a->len < b->len) {
        if (!bint_alloc(a, b->len))
            return NULL;
    }
    else {
        memset(a->digits + b->len, 0, (a->len - b->len) * sizeof(word));
    }

    memcpy(a->digits, b->digits, b->len * sizeof(word));
    a->len = b->len;
    a->sgn = b->sgn;
    return a;
}

bint *bint_add(bint *result, bint *a, bint *b) {
    return _bint_add_digits(result, a, b, b->sgn^a->sgn);
}

bint *bint_sub(bint *result, bint *a, bint *b) {
    return _bint_add_digits(result, a, b, !(b->sgn^a->sgn));
}

bint *bint_add_to(bint *a, bint *b) {
    return bint_add(a, a, b);
}

bint *bint_sub_to(bint *a, bint *b) {
    return bint_sub(a, a, b);
}
