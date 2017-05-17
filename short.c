#include <stdlib.h>
#include <memory.h>

#include "bint.h"
#include "atomics_inline.c"


bint *_bint_shrt_add    (bint *result, bint *a, word b, ulong flip);
bint *_bint_shrt_add_to (              bint *a, word b, ulong flip);


word braw_short_addmul(word *dst, word *src, ulong dst_len, ulong src_len, word mb) {
    word mc = 0, /* multiplication carry */
         ac = 0; /* addition carry */
    ASSERT(dst && src);
    ASSERT(dst_len >= src_len);
    dst_len -= src_len;
    do *dst = add_with_carry(*dst, mul_with_carry(*src++,mb,&mc), &ac), dst++;
     while (--src_len);
    if (dst_len && mc) {
        *dst = add_with_carry(*dst, mc, &ac);
        mc = 0, dst++, dst_len--;
    }
    while (dst_len-- && ac)
        *dst = add_with_carry(*dst, 0, &ac);
    return mc + ac;
}

word braw_short_submul(word *dst, word *src, ulong dst_len, ulong src_len, word mb) {
    word mc = 0, /* multiplication carry */
         sb = 0; /* subtraction borrow */
    ASSERT(dst && src);
    ASSERT(dst_len >= src_len);
    dst_len -= src_len;
    do *dst = sub_with_borrow(*dst, mul_with_carry(*src++,mb,&mc), &sb), dst++;
     while (--src_len);
    if (dst_len && mc) {
        *dst = sub_with_borrow(*dst, mc, &sb);
        dst++, dst_len--;
    }
    for (mc = 0; dst_len && sb; dst++,dst_len--)
        *dst = sub_with_borrow(*dst, 0, &sb);
    return mc + sb;
}


word braw_short_add(word *dst, word *src, ulong src_len, word b) {
    ASSERT(dst && src);
    if (!src_len--) return b;
    *dst = b + *src++;
    if (*dst++ < b) {
        while (src_len--) {
            if (*src == WORDMAX) {
                *dst++ = 0;
                 src++;
            } else {
                *dst++ = 1+*src++;
                if (src_len > 1)
                    memcpy(dst,src,(src_len-1)*sizeof(word));
                return 0;
            }
        }
        return 1;
    } else {
        memcpy(dst,src,src_len*sizeof(word));
        return 0;
    }
}
word braw_short_add_to(word *src, ulong src_len, word b) {
    ASSERT(src);
    if (!src_len--) return b;
    *src += b;
    if (*src++ < b) {
        while (src_len--) if (++(*src++)) return 0;
        return 1;
    } else return 0;
}


word braw_short_sub(word *dst, word *src, ulong src_len, word b) {
    word loword_overflow;
    ASSERT(dst && src);
    if (!src_len--) return 0-b;
    loword_overflow = b - *src;
    *dst++ = *src++ - b;
    if (loword_overflow < b) {
        while (src_len--) {
            if (!*src) {
                *dst++ = WORDMAX;
                 src++;
            } else {
                *dst++ = *src++ - 1;
                if (src_len > 1)
                    memcpy(dst, src, (src_len-1)*sizeof(word));
                return 0;
            }
        }
        return loword_overflow;
    } else {
        memcpy(dst, src, src_len*sizeof(word));
        return 0;
    }
}

word braw_short_sub_to(word *src, ulong src_len, word b) {
    word loword_overflow;
    ASSERT(src);
    if (!src_len--) return 0-b;
    loword_overflow = b - *src;
    *src++ -= b;
    if (loword_overflow < b) {
        while (src_len--) if ((*src++)--) return 0;
        return loword_overflow;
    } else return 0;
}


word braw_inc(word *dst, ulong len) {
    ASSERT(dst);
    while (len--) if (++(*dst++)) return 0;
    return 1;
}

word braw_dec(word *dst, ulong len) {
    ASSERT(dst);
    while (len--) if ((*dst++)--) return 0;
    return 1;
}

void braw_neg(word *dst, ulong len) {
    ASSERT(dst);
    while (len--) {
        *dst = (~*dst) + 1;
        if(*dst++) break;
    }
    while (len--) {
        *dst = ~*dst;
        dst++;
    }
}



word braw_short_mul(word *dst, word *src, ulong src_len, word b) {
    word carry = 0;
    ASSERT(dst && src);
    do *dst++ = mul_with_carry(*src++,b,&carry); while (--src_len);
    return carry;
}

INLINE word braw_short_mul_to(word *src, ulong src_len, word b) {
    return braw_short_mul(src,src,src_len,b);
}


word braw_short_div(word *dst, word *src, ulong src_len, word b) {
    word remainder = 0;
    ASSERT(dst && src);

    src += src_len - 1;
    dst += src_len - 1;

    while (src_len--)
        *dst-- = div_with_remainder(*src--,b,&remainder);

    return remainder;
}

INLINE word braw_short_div_to(word *src, ulong src_len, word b) {
    return braw_short_div(src,src,src_len,b);
}



bint *_bint_shrt_add(bint *result, bint *a, word b, ulong flip) {
    if (a->sgn == flip) {
        if (!bint_alloc(result,a->len+1))
            return NULL;
        else result->len = a->len;
        if (result->digits[a->len] = braw_short_add(
                result->digits,a->digits,a->len,b))
            result->len++;
        result->sgn = a->sgn;
    } else {
        if (!bint_alloc(result,a->len))
            return NULL;
        else {
            word borrow;
            if (borrow = braw_short_sub(
                    result->digits, a->digits, a->len, b)) {
                result->digits[0] = borrow;
                result->sgn = a->sgn ^ 1;
                result->len = 1;
            } else {
                result->sgn = a->sgn;
                result->len = a->len - 1;
                if (!result->digits[result->len])
                    result->len--;
            }
        }
    }
    return result;
}

bint *bint_short_add(bint *r, bint *a, word b) { return _bint_shrt_add(r,a,b,0); }
bint *bint_short_sub(bint *r, bint *a, word b) { return _bint_shrt_add(r,a,b,1); }


bint *_bint_shrt_add_to(bint *a, word b, ulong flip) {
    if (a->sgn == flip) {
        if (!bint_alloc(a,a->len+1)) return NULL;
        if (a->digits[a->len] = braw_short_add_to(a->digits,a->len,b))
            a->len++;
    } else {
        word borrow;
        if (borrow = braw_short_sub_to(a->digits, a->len, b)) {
            a->digits[0] = borrow;
            a->sgn ^= 1;
            a->len  = 1;
        } else if (!a->digits[a->len-1]) a->len--;
    }
    return a;
}

bint *bint_short_add_to(bint *a, word b) { return _bint_shrt_add_to(a,b,0); }
bint *bint_short_sub_to(bint *a, word b) { return _bint_shrt_add_to(a,b,1); }



bint *bint_short_mulu(bint *result, bint *a, word b) {
    ASSERT(result && a);
    result->sgn = a->sgn;
    if (!bint_alloc(result,a->len + 1))
        return NULL;
    else result->len = a->len;
    if (result->digits[a->len] = braw_short_mul(
            result->digits,a->digits,a->len,(word)b))
        result->len++;
    return result;
}
bint *bint_short_mulu_to(bint *a, word b) {
    ASSERT(a);
    if (!bint_alloc(a,a->len+1)) return NULL;
    else if (a->digits[a->len] = braw_short_mul_to(a->digits,a->len,b))
        a->len++;
    return a;
}
bint *bint_short_mul(bint *result, bint *a, sint b) {
    ASSERT(result && a);
    if (b < 0) b = -b, result->sgn = a->sgn ^ 1;
    else result->sgn = a->sgn;
    return bint_short_mulu(result,a,(word)b);
}
bint *bint_short_mul_to(bint *a, sint b) {
    ASSERT(a);
    if (b < 0) b = -b, a->sgn ^= 1;
    return bint_short_mulu_to(a,(word)b);
}



bint *bint_short_divu(bint *result, bint *a, word b, word *remainder) {
    word r;
    ASSERT(result && a);   
    if (!bint_alloc(result,a->len)) 
        return NULL;
    r = braw_short_div(result->digits, a->digits, a->len, b);
	result->len = a->len;
	result->sgn = a->sgn;
	bint_fit(result);
    if (*remainder) 
        *remainder = r;
    return result;
}
bint *bint_short_divu_to(bint *a, word b, word *remainder) {
    word r;
    ASSERT(a);
    r = braw_short_div_to(a->digits,a->len,b);
	bint_fit(a);
    if (remainder) *remainder = r;
    return a;
}


bint *bint_short_div(bint *result, bint *a, sint b, sint *remainder) {
    word rmdr;
    if (b > 0) {
        if (bint_short_divu(result,a,(word)b,&rmdr)) {
            if (remainder) *remainder = (sint) rmdr;
            return result;
        } else return NULL;
    } else {
        bint aI = *a;
        aI.sgn ^= 1;
        if (bint_short_divu(result,&aI,(word)(-b),&rmdr)) {
            if (remainder) *remainder = -((sint) rmdr);
            return result;
        } else return NULL;
    }
}

bint *bint_short_div_to(bint *a, sint b, sint *remainder) {
    word rmdr;
    if (b > 0) {
        if (bint_short_divu_to(a,(word)b,&rmdr)) {
            if (remainder) *remainder = (sint) rmdr;
            return a;
        } else return NULL;
    } else {
        a->sgn ^= 1;
        if (bint_short_divu_to(a,(word)(-b),&rmdr)) {
            if (remainder) *remainder = -((sint)rmdr);
            return a;
        } else return NULL;
    }
}


bint* bint_short_assign(bint *a, sint b) {
    a->sgn = b<0?1:0;
    memset(a->digits,0,a->len*sizeof(word));
    a->digits[0] = (word) ABS(b);
    a->len = b ? 1 : 0;
    return a;
}


bint *bint1(bint *a) {
	return bint_short_assign(a, 1);
}

bint *bint0(bint *a) {
	return bint_short_assign(a, 0);
}
