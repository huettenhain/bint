#ifndef _EYE_HPP
#define _EYE_HPP

#include <float.h>
#include "macros.h"

typedef unsigned long  ulong;
typedef unsigned short ushort;

#ifndef __cplusplus
 typedef int bool;
 const static int true  = 1;
 const static int false = 0;
#endif

#if defined(BITS64)
  typedef    int64_t    sint;
  typedef   uint64_t    word;
# define    WORDSIZE    0x40
# define     WORDMAX    0xFFFFFFFFFFFFFFFF
#elif defined(BITS32)
  typedef    int32_t    sint;
  typedef   uint32_t    word;
# define    WORDSIZE    0x20
# define     WORDMAX    0xFFFFFFFF
#else
# include  <limits.h>
  typedef signed long   sint;
  typedef       ulong   word;
# define     WORDSIZE   (sizeof(word)/sizeof(char)*CHAR_BIT)
# define      WORDMAX   UINT_MAX
#endif


#define WONE   ((word)1)
#define HWSIZE (WORDSIZE >> 1)
#define HWMASK ((WONE<<HWSIZE)-1)
#define LO(_w) ((_w) & HWMASK)
#define HI(_w) ((_w) >>HWSIZE)

#define ABS(_b) (((_b)>0)?(_b):(-(_b)))

#define BINT_MSB_MASK (WONE << (WORDSIZE-1))
#define BINT_INITSIZE (0x200/WORDSIZE)


typedef struct {
    ulong bits_per_double;
    ulong length_of_vector;
} fmul_context;

typedef struct bmc {  /* bigint multiplication context */
    word  *buffer;    /* allocated buffer */
    ulong  size;      /* size of allocated buffer */

 /* any integer with less than k_cut many words is multiplied classically,
    anything above it will use the divide-and-conquer method. */
    ulong  k_cut;
#   define BINT_KCUT 15

 /* For each bigint a, let B=(1<<i) be the smallest power of two such that 
    the length of a is less than B. Then, if the length of a is at most f_dst[i]
    words less than B, we deploy the FFT multiplication. */
    sint   f_dst[WORDSIZE - 1];
} bmc;




typedef struct bint {
    bmc    *context;  /* thread-specific context (multiplication buffer) */
    ulong   len;      /* highest index of a nonzero digit */
    ulong   mem;      /* number of words allocated */
    word   *digits;   /* array of digits of the big integer */
    ulong   sgn;      /* sign bit. 1:negative, 0:positive */
} bint;

bmc  *bmc_create        ();
void  bmc_destroy       (bmc *context);
bool  bmc_resize        (bmc *context, ulong bytes);

void  bmc_defaults      (bmc *context); /* load default values */
void  bmc_fft           (bmc *context); /* always use FFT method */

bint *bint_init         (bmc *context, sint init);
bint *bint_init_alloc   (bmc *context, sint init, ulong size);

/* reallocate the memory of a. returns NULL if not enough memory is available. */
bint *bint_alloc        (bint *a, ulong size);
void  bint_free         (bint *a);

/* after a division, reduce a->len until it is correct */
void  bint_fit          (bint *a);
bint *bint_from_str     (bint *a, char *s);
char *bint_to_str       (bint *a, ulong base);


/*******************************************************************************

 BigInt Library General Notes

 Abstract Functions and Raw Functions:
 Every function of the big int library is prefixed with either one of bint_ or
 braw_. We will refer to them as abstract and raw functions, respectively. The 
 differences are the following:

 - Abstract functions operate on bint structures and allocate memory when it
   is necessary to grow the digits array. Abstract function names are sometimes
   suffixed with a single character 'u' to indicate that the function operates
   on the unsigned part of the bigint. This is sometimes useful if you do not
   require signed numbers, and you want to avoid the additional checks that 
   become necessary when using signed numbers.

 - Raw functions operate on word-arrays which contain the digits representing
   an unsigned big int. We refer to these operands as raw big ints. They are 
   given by a pointer to the array of digits and an ulong value containing the
   size of this array (in words). 
   In the documentation of a raw function, a number is referred to by the 
   syntax (ptr,len) where
     - ptr is the name of a pointer to the digits of the big int
     - len the name of an ulong variable containing the number of digits at ptr
   NOTE: Using this syntax in documentation, we implicitly assume that there is
         enough space at ptr to store len many words. 
   NOTE: Raw functions NEVER allocate memory.

********************************************************************************
 
 0. Base and Helper Functions

********************************************************************************

 0.1. Arithmetic Base Functions 

      These functions represent the basic arithmetic operations provided by
      the CPU itself. 

      Since the functions are inlined, they don't need to have a declaration
      here.

word div_with_remainder(word a, word b, word *r);
word add_with_carry(word a, word x, word *c);
word mul_with_carry(word a, word b, word *c);
word sub_with_borrow(word a, word x, word *b);


 0.2. Helper Functions
      A set of non-arithmetic functions that are frequently required.

      Implementation: bint.c

*******************************************************************************/


/* The comparison routine for raw bigints yields
     -1 ; if (a,n) is smaller than (b,m).
      0 ; if (a,n) is equal to (b,m).
      1 ; if (a,n) is greater than (b,m). */
int braw_compare (word *a, word *b, ulong a_len, ulong b_len);
int bint_compare (bint *a, bint *b);
bool bint_is_zero(bint *a);
bool bint_is_sint(bint *a);
sint bint_to_sint(bint *a);

/* return the logarithm of (a,n) to the base B=(1<<WORDSIZE), rounded down */
ulong braw_logB(word *a, ulong n); 

/* return the logarithm of (a,n) to the base 2, rounded down. */
ulong braw_log2(word *a, ulong n);

ulong   bint_logB(bint *a);     /* log_B(a), rounded down. B=(1<<WORDSIZE) */
ulong   bint_log2(bint *a);     /* log_2(a), rounded down. */

ulong   word_log2(word a);      /* log_2(a), rounded down. */
ulong   word_log2_ceil(word a); /* the smallest n such that a <= 2^n */
word    word_next_pow2(word a); /* 2^word_log2_ceil(a) */
bool    word_is_pow2(word a);   /* nonzero <=> a is a power of 2 */

/*******************************************************************************

 1. Arithmetic Operations

    The following functions implement the basic binary arithmetic operations
    on big ints:
     
      - addition
      - subtraction
      - multiplication
      - integer division with remainder

    Abstractly spoken, these functions always have three paramters: Two 
    source parameters (which are to be composed) and one result parameter
    (which is supposed to store the result of the composition). 
    
    Functions with a _to suffix (from hereon called TO-functions) always
    work exactly as their equivalent without suffix, except that the result
    operand is the same as the first source operand. 
    
    TO-Functions should be applied whenever possible for optimal performance.

********************************************************************************

 1.1. Short Arithmetic Operations.

      Here, the first source operand is a big integer while the second source 
      operand is only one machine word in size. 

      Implementation: short.c

********************************************************************************

 1.1.1 Raw short arithmetic functions: Perform an arithmetic composition of an
       unsigned raw big int and a short operand.

*******************************************************************************/

/* Add b to (src,src_len) and store the result in (dst,src_len). If src_len is 
   zero, the function returns b. Otherwise, the carrybit is returned. */
word braw_short_add    (word *dst, word *src, ulong src_len, word b);
word braw_short_add_to (           word *src, ulong src_len, word b);

/* Subtract b from (src,src_len) and store the result in (dst,src_len). If an
   overflow occurs, the function returns -dst[0]. If src_len is 0, the function
   returns 0-b. In all other cases, the return value is zero. */
word braw_short_sub    (word *dst, word *src, ulong src_len, word b);
word braw_short_sub_to (           word *src, ulong src_len, word b);

/* Multiply (src,src_len) with the word b and store the result in (dst,src_len).
   The function returns the carryword that occurs. */
word braw_short_mul    (word *dst, word *src, ulong src_len, word b);
word braw_short_mul_to (           word *src, ulong src_len, word b);

/* Performs a nonnegative integer division of (src,src_len) by the word b
   and stores the quotient in (dst,src_len). If b is zero, a division by zero
   will occur and behaviour is undefined. The function returns the word-sized
   remainder of the division. */
word braw_short_div    (word *dst, word *src, ulong src_len, word b);
word braw_short_div_to (           word *src, ulong src_len, word b);

/* The following methods first multiply the src operand with the word mb and
   then add/subtract it to/from the destination operand. The destination 
   operand is assumed to have at least (src_len+1) words allocated. 
   The carry/borrow word of the operation is returned. */
word braw_short_addmul(word *dst, word *src, ulong dst_len, ulong src_len,
                       word mb);
word braw_short_submul(word *dst, word *src, ulong dst_len, ulong src_len, 
                       word mb);

/* methods to increase/decrease a raw big integer. */
word braw_inc(word *dst, ulong len);
word braw_dec(word *dst, ulong len);


/*******************************************************************************

 1.1.2 Abstract short arithmetic functions: Perform an arithmetic composition
       of a signed big int with a second, signed, one-word operand. 

*******************************************************************************/

/* Addition, subtraction and multiplication with short operands work in the 
   obvious way. */
bint *bint_short_add    (bint *result, bint *a, word b);
bint *bint_short_add_to (              bint *a, word b);
bint *bint_short_sub    (bint *result, bint *a, word b);
bint *bint_short_sub_to (              bint *a, word b);

/* Multiplication with signed one-word operands */
bint *bint_short_mul    (bint *result, bint *a, sint b);
bint *bint_short_mul_to (              bint *a, sint b);
/* Multiplication with unsigned one-word operands */
bint *bint_short_mulu   (bint *result, bint *a, word b);
bint *bint_short_mulu_to(              bint *a, word b);

/* Division by signed one-word operands. The short division function takes an
   optional parameter remainder in which the remainder of the integer division
   is stored. Recall how integer division is defined:

   Let a and b be integers with b not equal to zero. Then, the inger division
   of a by b is given by the integers q and r such that 
     a = bq + r  and |r| < |b|.
   If q is minimal with this property, we write div(a,b) := (q,r). Example:

     div(  23,  4 ) = (  5,  3 )                                          (SDIV)
     div( -23, -4 ) = (  5, -3 )

   and

     div( -23,  4 ) = ( -6,  1 )
     div(  23, -4 ) = ( -6, -1 )
     
   Another technical comment - the remainder parameter may be NULL. In this
   case, it is ignored. */
bint *bint_short_div    (bint *result, bint *a, sint b, sint *remainder);
bint *bint_short_div_to (              bint *a, sint b, sint *remainder);

/* Division by unsigned one-word operands, ignores the sign completely. */
bint *bint_short_divu   (bint *result, bint *a, word b, word *remainder);
bint *bint_short_divu_to(              bint *a, word b, word *remainder);


/*******************************************************************************

 1.2 Long Arithmetic Operations.

     In long arithmetic operations, both the first and the second operand
     are big integers. 

********************************************************************************

 1.2.1 Raw Long Additive Operations
       Implementation: long.additive.c
       
*******************************************************************************/

/* Add the bigint b to the bigint a and store the result in dst. It is assumed
   that b_len is less or equal to a_len and that at dst, at least a_len cells
   of memory are available. The carry of the operation is returned. */
word braw_shorter_add   (word *dst, word *a, word *b, ulong a_len, ulong b_len);
word braw_shorter_add_to(           word *a, word *b, ulong a_len, ulong b_len);

/* Subtract b from a and store the result in dst. It is assumed that b_len is
   less or equal to a_len and that at dst, at least a_len cells of memory are
   available. The borrow of the operation is returned. */
word braw_shorter_sub   (word *dst, word *a, word *b, ulong a_len, ulong b_len);
word braw_shorter_sub_to(           word *a, word *b, ulong a_len, ulong b_len);

/* Add / Subtract two raw unsigned bigints of the same length. The carry /
   borrow of the operation is returned. */
word braw_add    (word *dst, word *a, word *b, ulong len);
word braw_add_to (           word *a, word *b, ulong len);

word braw_sub    (word *dst, word *a, word *b, ulong len);
word braw_sub_to (           word *a, word *b, ulong len);

void braw_neg    (word *dst, ulong len);


/*******************************************************************************

 1.2.2 Raw Long Multiplicative Operations
       Implementation: long.multiplicative.c
       
*******************************************************************************/


/* Grade school multiplication.
  
   Input
     (a,a_len) raw bigint
     (b,b_len) raw bigint
     (dst,a_len+b_len) buffer
   Output
     (dst,a_len+b_len) = a*b
*/   
void braw_gmul(
    word *dst, 
    word *a, 
    word *b, 
    ulong a_len, 
    ulong b_len
);


/* Karatsuba-Offman multiplication.
  
   Input
     (a,a_len) raw bigint
     (b,b_len) raw bigint
     (dst,KMUL_BUFFERSIZE(MAX(a_len,b_len))) buffer
   Output
     (dst,a_len+b_len) = a*b
*/
void braw_kmul(
        word *dst, 
        word *a, 
        word *b, 
        ulong a_len, 
        ulong b_len, 
        ulong cut
);
#define KMUL_BUFFERSIZE(n) (  \
    sizeof(word) * (16 * (word_log2(n)+1) +  6 * n  ))


void  braw_fmul( fmul_context *ctx, 
                 void* bfr, 
                 word *dst, 
                 word *a, 
                 word *b, 
                 ulong a_len, 
                 ulong b_len);

void  braw_fmul_getcontext( fmul_context *ctx, ulong a_len, ulong b_len );
ulong braw_fmul_buffersize( fmul_context *ctx );


/* Grade school division with scaling.  ***********
   Input
     (a,a_len) raw bigint
     (b,b_len) raw bigint
     (dst,GDIV_BUFFERSIZE(a_len,b_len)) buffer
   Output
     return k
     If k = 0, then a < b meaning a/b=0 and a%b=a.
     Else:
     (dst,k)        = a/b
     (dst+k, b_len) = a%b
     
***************************************************/
ulong braw_divmod(
    word *dst,
    word *a,
    word *b,
    ulong a_len,
    ulong b_len
);
#define GDIV_BUFFERSIZE(_n,_m) (      \
     ((2*(_n))-(_m)+2) * sizeof(word) \
)

/*******************************************************************************

 1.2.3 Abstract Long Unsigned Arithmetic Operations: These are all the abstract
       arithmetic operations which ignore the sign.

*******************************************************************************/




/*******************************************************************************

 1.2.4 Abstract Long Signed Arithmetic Operations: This is as straight forward
       as it gets, here we have all arithmetic operations that take only big
       integers as operands.

       Implementation: long.additive.c
                       long.multiplicative.c
       
*******************************************************************************/

bint *bint_assign      (bint *a, bint *b);
bint *bint_short_assign(bint *a, sint  b);

bint *bint1(bint *a);
bint *bint0(bint *a);

bint *bint_add    (bint *result, bint *a, bint *b);
bint *bint_add_to (              bint *a, bint *b);

bint *bint_sub    (bint *result, bint *a, bint *b);
bint *bint_sub_to (              bint *a, bint *b);

bint *bint_mul    (bint *result, bint *a, bint *b);
bint *bint_mul_to (              bint *a, bint *b);


/* The division works as in (SDIV) and is defined in terms of the routine
   bint_divmod. It can return one of the following values: */
typedef enum  {
    BDIV_OK      = 0, /* division successful */
    BDIV_ZERODIV = 1, /* attempt of a zero division */
    BDIV_MEMFAIL = 2, /* a memory allocation failed */
} bint_divmod_result;

bint_divmod_result bint_divmod(bint *result, bint *a, bint *b, bint *remainder);

/* This function can be used to calculate the quotient and the remainder 
   at once. In terms of bint_divmod, the functions bint_div and bint_mod are 
   defined, which allow you to return only the quotient or the remainder. */

bint *bint_div    (bint *result, bint *a, bint *b);
bint *bint_div_to (              bint *a, bint *b);

bint *bint_mod    (bint *result, bint *a, bint *b);
bint *bint_mod_to (              bint *a, bint *b);


/*******************************************************************************

 2. Shifting Functions

    All necessary binary shifting operators have been implemented. Shifting 
    by one single bit has been implemented separately. This is mostly done
    for convenience, the performance advantage over passing 1 as the places
    parameter to the more general functions should be minimal. On the other 
    hand, this has not been tested.

    Remark: The bint_lshift routine will allocate space for the shift and fail
    if no memory can be allocated. Shifting to the right will always succeed.

    Implementation: shifts.c

*******************************************************************************/

/* return number of words nulled */
ulong braw_lshift  (word *a, ulong size, ulong places);
ulong braw_rshift  (word *a, ulong size, ulong places);
void  braw_rshift1 (word *a, ulong size);
void  braw_lshift1 (word *a, ulong size);

bint *bint_lshift  (bint *a, ulong places);
bint *bint_rshift  (bint *a, ulong places);
bint *bint_lshift1 (bint *a);
bint *bint_rshift1 (bint *a);


#endif /***********************************************************[ eof ]*****/

