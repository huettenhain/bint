#include <string.h>
#include <math.h>

#include "bint.h"
#include "fft.h"
#include "base.c"


void braw_gmul(word *dst, word *a, word *b, ulong a_len, ulong b_len)
{
    ulong k;
    ASSERT( dst && a && b  );
    ASSERT( dst != a && dst != b );
    ASSERT( b_len <= a_len );
    dst[a_len] = braw_short_mul(dst,a,a_len,b[0]);
    for (k=1; k<b_len; k++)
        dst[a_len+k] = braw_short_addmul(dst+k,a,a_len,a_len,b[k]);
}



void braw_kmul(word *dst, word *a, word *b, ulong n, ulong m, ulong cut) {
 /* n = length of a,
    m = length of b. */
    ASSERT( dst && a && b );
    ASSERT( dst != a && dst != b);
    ASSERT( m <= n );
    if (n <= cut) {  
        braw_gmul(dst,a,b,n,m);
    } else {
        ulong k = n/2;
        ulong l = n-k;
        if (m <= k) { /* easy split, b is smaller than half the size of a */
            braw_kmul(dst,     a,   b, k, m, cut);  /* dst[0] - dst[k+m-1] = a0*b */
            braw_kmul(dst+n+m, a+k, b, l, m, cut);  /* bfr[0] - bfr[l+m-1] = a1*b */
            braw_shorter_add(dst+k, dst+n+m, dst+k, l+m, m);
        } else {
            ulong j = m-k, t;
            braw_kmul(dst,     a,   b,   k, k, cut);  
            braw_kmul(dst+2*k, a+k, b+k, l, j, cut);   /* dst[0] - dst[m+n-1] = T2|T0 */
            memcpy(dst+m+n, dst, 2*k*sizeof(word));  /* backup of T0 */
            braw_shorter_sub_to(dst+k, dst+2*k, m+l, l+j); /* center minus T2 */
            braw_shorter_sub_to(dst+k, dst+m+n, m+l, 2*k); /* center minus T0 */    
            /* now we calculate T1 = (a0+a1)(b0+b1) and add it to center. */
            dst += n+m;
            dst[l] = braw_shorter_add(dst, a+k, a, l, k);
            if (j>=k)  
                dst[l+1+j] = braw_shorter_add(dst+l+1, b+k, b, t=j, k);
            else     
                dst[l+1+k] = braw_shorter_add(dst+l+1, b, b+k, t=k, j);
            braw_kmul(dst+l+t+2, dst, dst+l+1, l+1, t+1, cut);
            braw_shorter_add_to(dst-l-m, dst+l+t+2, m+l, MIN(l+t+2,m+l));
        }
    }
}


INLINE void _feed_word_to_double( double *d, word *s, ulong N, ulong n, ulong bit_stuffing) {
    ulong bits_left,i,j;
    word tmp;
    for (bits_left=WORDSIZE,i=0,j=1,tmp=*s; i<N; i++) {
        if (bits_left < bit_stuffing) {
            ulong bits_needed = bit_stuffing - bits_left;
            d[i] = (double) tmp;
            if (j < n) {
                tmp = s[j++];
                d[i] += ( tmp & ((1<<bits_needed)-1) ) << bits_left;
                tmp >>= bits_needed;
                bits_left = WORDSIZE - bits_needed;
            }
        } else {
            d[i] = (double) (tmp & ((1<<bit_stuffing)-1));
            tmp >>= bit_stuffing;
            bits_left -= bit_stuffing;
        } 
    }
}


INLINE void _feed_double_to_word( word *d, double *s, ulong n, ulong bit_stuffing) {
    ulong bits_in_buffer = bit_stuffing,
          bits_left = WORDSIZE, i=0, j=1;
    word tmp=(word)s[0];
    for (*d=0;;) {
        word bits_occupied = WORDSIZE - bits_left;
        if (bits_left < bits_in_buffer) {
            d[i++] |= ( tmp & ((1<<bits_left)-1) ) << bits_occupied;
            if (i == n) return;
            d[i] = 0;
            tmp >>= bits_left;
            bits_in_buffer -= bits_left;
            bits_left = WORDSIZE;
        } else {
            d[i] |= tmp << bits_occupied;
            tmp = (word) s[j++];
            bits_left -= bits_in_buffer;
            bits_in_buffer = bit_stuffing;
        }
    }
}


void braw_fmul_getcontext( fmul_context *ctx, ulong n, ulong m ) {
    ulong result_size;
#   define l ctx->bits_per_double
#   define k ctx->length_of_vector
    
    result_size = (n+m)*WORDSIZE;
    l = LDBL_MANT_DIG/2-1;
    k = word_log2_ceil(DIVUP(result_size, l));

    while (LDBL_MANT_DIG < 2*l + k + word_log2(k)*(word_log2(k)+1)) {
        l = l - 1;
        k = word_log2_ceil(DIVUP(result_size, l));
    }
#   undef l
#   undef k
}

ulong braw_fmul_buffersize( fmul_context *ctx ) {
    return TWO_TO_THE(ctx->length_of_vector) * sizeof(double) * 2;
}

void braw_fmul(fmul_context *ctx, void *bfr, word *dst, word *a, word *b, ulong n, ulong m) {

    double *A,*B,T,C;
    ulong N,M,i,size;

    ASSERT( bfr && a && b );
    ASSERT( bfr != a && bfr != b);

    N = DIVUP( n * WORDSIZE, ctx->bits_per_double );
    M = DIVUP( m * WORDSIZE, ctx->bits_per_double );

    size = 1 << ctx->length_of_vector;
    
    A = (double*) bfr;
    B = A+size;

    _feed_word_to_double(A,a,N,n,ctx->bits_per_double);
    for (i=N;i<size;i++) A[i] = 0.0;

    _feed_word_to_double(B,b,M,m,ctx->bits_per_double);
    for (i=M;i<size;i++) B[i] = 0.0;
    
    realfft((COMPLEX*)A,size/2,-1);
    realfft((COMPLEX*)B,size/2,-1);

    A[0] *= B[0]; A[1] *= B[1];
    for (i=2;i<size;i+=2) {
        A[i]   = (T=A[i])*B[i]-A[i+1]*B[i+1];
        A[i+1] =  T*B[i+1] + A[i+1]*B[i];
    }

    realfft((COMPLEX*)A,size/2,1);

#   define CARRYFIX ((double)(1<<ctx->bits_per_double))

    for (C=0.0,i=0; i<size; i++) {
        C = floor( (T = A[i]*size+C+.5) / CARRYFIX );
        A[i] = T - C*CARRYFIX;
    }

#   undef CARRYFIX

    _feed_double_to_word( dst, A, n+m, ctx->bits_per_double );
}


ulong braw_divmod(word *dst, word *a, word *b, ulong n, ulong m) {

    ulong d,k,j;
    word carry = 0, *r;
    
    ASSERT( dst && a && b );
    ASSERT( dst != a && dst != b );
    ASSERT( n >= m );
    
    d = WORDSIZE - word_log2(b[m-1]) - 1;
    k = n-m+1;

    memcpy( (r = dst+k), a, n * sizeof(word));
    r[n] = 0;

    braw_lshift(r, n+1, d); /* scale the divident */
    braw_lshift(b,   m, d); /* scale b */

    for (j = n-m; j < n; j--) {
        dst[j] = (r[j+m] >= b[m-1]) ? WORDMAX :
            div_without_remainder(r[j+m],r[j+m-1],b[m-1]); /* guess the next digit d */  
        if (braw_short_submul(r+j,b,m+1,m,dst[j]))
            do dst[j]--;
                while (!braw_shorter_add_to(r+j,b,m+1,m)); /* digit was too large */
    }

    braw_rshift(r, n+1, d);  /* unscale the remainder */
    braw_rshift(b,   m, d);  /* unscale b */

    return k;
}



bint *bint_mul(bint *result, bint *a, bint *b) {   
    sint l2;
    ASSERT(result && a && b);
       
	if (a->len < b->len)
		return bint_mul(result, b, a);
	else if (bint_is_zero(b))
		return bint0(result);
   
    if (!bint_alloc(result, a->len + b->len))
        return NULL;

    l2 = word_log2_ceil(a->len);
	ASSERT(l2 < WORDSIZE - 1); /* this is quite reasonable. */

    if ((TWO_TO_THE(l2) - (sint)a->len) < result->context->f_dst[l2]) {
        fmul_context ctx;
        braw_fmul_getcontext( &ctx, a->len, b->len );

        if (bmc_resize(result->context, braw_fmul_buffersize( &ctx ))) {        
            braw_fmul(
                &ctx,
                result->context->buffer,
                result->digits,
                a->digits,
                b->digits,
                a->len,
                b->len
            );
            result->len = a->len + b->len;
            goto _bint_mul_done;
        }
    } 

    
    if (bmc_resize(result->context, KMUL_BUFFERSIZE(a->len))) {
        braw_kmul(
            result->context->buffer,
            a->digits,
            b->digits,
            a->len,
            b->len,
            result->context->k_cut
        );
        memcpy( result->digits, result->context->buffer, 
            (result->len = a->len + b->len) * sizeof(word));
    } else {
        return NULL;
    }

_bint_mul_done:
    result->len = bint_logB(result)+1;
    result->sgn = a->sgn ^ b->sgn;
    return result;
}


bint *bint_mul_to(bint *a, bint *b) {
    return bint_mul(a,a,b);
}


bint_divmod_result bint_divmod(bint *result, bint *a, bint *b, bint *remainder) {  
    bmc *context;
	ulong new_sgn;

    ASSERT(a && b);
    ASSERT(result || remainder);

    if (bint_is_zero(b))
        return BDIV_ZERODIV;
       
    switch( braw_compare( a->digits, b->digits, a->len,  b->len ) ) {
    case  0: 
        if (remainder)
            bint_short_assign( remainder, 0 );
        if (result) 
            if (a->sgn ^ b->sgn) 
                bint_short_assign( result, -1 );
            else
                bint_short_assign( result,  1 );
        return BDIV_OK;
    case -1:
        if (a->sgn ^ b->sgn) {
            if (remainder && !bint_add(remainder,b,a))
                return BDIV_MEMFAIL;
            if (result)
                bint_short_assign(result,-1);
        } else {
            if (remainder && !bint_assign(remainder,a))
                return BDIV_MEMFAIL;
            if (result)
                bint_short_assign(result,0);
        }
        return BDIV_OK;
    case  1: 
        context = result ? result->context : remainder->context;
   
        if (bmc_resize(context, GDIV_BUFFERSIZE( a->len, b->len ))) {   
            ulong r,l;

            r = braw_divmod(
                context->buffer,
                a->digits,
                b->digits,
                a->len,
                b->len);
    
            if (result) {
                l = braw_logB( context->buffer, r ) + 1;
                if (!bint_alloc(result,l)) return BDIV_MEMFAIL;
                else memcpy( result->digits, context->buffer, 
                    (result->len = l) * sizeof(word) );
            }

            if (remainder) {
                l = braw_logB( context->buffer + r, b->len ) + 1;
                if (!bint_alloc(remainder,l)) return BDIV_MEMFAIL;
                else memcpy( remainder->digits, context->buffer + r,
                    (remainder->len = l)*sizeof(word) );
                remainder->sgn = b->sgn;
            }

            if ( (new_sgn = a->sgn ^ b->sgn) ) {
                if (result) {
					if (braw_short_add_to(result->digits, result->len, 1)) {
						/* if a carry occurs here, it can only be a carry of 1 */
						if (!bint_alloc(result, result->len + 1))
							return BDIV_MEMFAIL;
						else result->digits[result->len++] = 1;
					}
                }
                if (remainder) {
                    braw_shorter_sub(
                        remainder->digits,
                        b->digits,
                        remainder->digits, 
                        b->len,
                        remainder->len);
                    remainder->len = bint_logB( remainder ) + 1;
                }
            }

			if (result)
				result->sgn = new_sgn;

            return BDIV_OK;
        } else {
            return BDIV_MEMFAIL;
        }

    NODEFAULT; }
}

bint *bint_div(bint *result, bint *a, bint *b) {
	ASSERT(!bint_is_zero(b));
    return (bint_divmod(result,a,b,NULL) == BDIV_OK) ? result : NULL;
}

bint *bint_div_to (bint *a, bint *b) {
	return bint_div(a, a, b);
}

bint *bint_mod(bint *result, bint *a, bint *b) {
	ASSERT(!bint_is_zero(b));
    return (bint_divmod(NULL,a,b,result) == BDIV_OK) ? result : NULL;
}

bint *bint_mod_to (bint *a, bint *b) {
	return bint_mod(a, a, b);
}


