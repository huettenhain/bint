#include "bint.h"

/* Calculate a + x + c where c is either 0 or 1 (carrybit).
   Returns the result and sets c to the carrybit of the
   addition. Subtraction with borrow is done similarly. */
INLINE word add_with_carry(word a, word x, word *c) {
#if defined(MVC_X86_64)
	word re;
	*c = _addcarry_u64(*c, a, x, &re);
	return re;
#else
    if (*c && !++a) 
        return x;
    *c = (x>WORDMAX-a) ? 1 : 0;
    return a+x;
#endif
}

INLINE word sub_with_borrow(word a, word x, word *b) {
#if defined(MVC_X86_64)
	word re;
	*b = _subborrow_u64(*b, a, x, &re);
	return re;
#else
    if (*b && !(a--))
        goto _carry_ok;
    *b = (a<x) ? 1 : 0;
_carry_ok:
    return a-x;
#endif
}

/* Calculate (a*b)+(*c). Return the lower word of the calculation
   and set (*c) to the upper word. */
INLINE word mul_with_carry(word a, word b, word *c) {
#if defined(MVC_X86_32)
    __asm {
        mov ecx, c
        mov eax, a
        mul b
        add eax, [ecx]
        adc edx, 0
        mov [ecx], edx
    }
#elif defined(MVC_X86_64)
	word hi, lo;
	a = _umul128(a, b, &hi);
	hi += _addcarry_u64(0, a, *c, &lo);
	*c = hi;
	return lo;
#elif defined(GCC_X86_32)
    word ret;
    asm( "mull %4"          "\n\t"
         "addl %0,%%eax"    "\n\t"
         "adcl $0,%%edx"    "\n\t"
         "movl %%edx,%0"    "\n\t" 
         "movl %%eax,%1"
    : "=g"(*c),"=g"(ret)
    :  "0"(*c), "a"(a),"g"(b)
    : "edx"
    );
    return ret;
#elif defined(GCC_X86_64)
    word ret;
    asm( "mulq %4"          "\n\t"
         "addq %0,%%rax"    "\n\t"
         "adcq $0,%%rdx"    "\n\t"
         "movq %%rdx,%0"    "\n\t" 
         "movq %%rax,%1"
    : "=g"(*c),"=g"(ret)
    :  "0"(*c), "a"(a),"g"(b)
    : "rdx"
    );
    return ret;
#else
    word al = LO(a), ah = HI(a),
         bl = LO(b), bh = HI(b);
    /* We calculate a*b by splitting a and b into one high and 
       one low half-word */

    a = al * bl; /* calculate preliminray lower word */
    b = ah * bh; /* calculate preliminary upper word */
    
    /* Calculate al * bh. We need to add it, shifted by HWSIZE
       bits, to the double-word [b][a]. To do so, we add the 
       upper half unshifted to [b]. The lower half is shifted
       and will be added to [a] later. */
    al *= bh; b += HI(al); al <<= HWSIZE;
    bl *= ah; b += HI(bl); bl <<= HWSIZE;

    if (*c > WORDMAX - a) b++; a += *c; /* add carry */

    if (al > WORDMAX - a) b++; a += al; /* add lower half of al*bh */
    if (bl > WORDMAX - a) b++; a += bl; /* add lower half of bl*ah */

    *c = b; return a;
#endif
}


#if 1 || !defined(MVC_X86_64)

/* Calculate [r][a] / [b]. Return the quotient and store
   the new remainder in [r]. */
INLINE word div_with_remainder(word a, word b, word *r) {
#if defined(MVC_X86_32)
    __asm {
        mov ebx,r
        mov eax,a
        mov edx,[ebx]
        div b
        mov [ebx],edx
    }
#elif defined(GCC_X86_32)
    word ret;
    asm( "divl %4" "\n\t" "movl %%eax,%1"
    : "=d"(*r),"=g"(ret)
    :  "0"(*r), "a"(a),"g"(b)
    );
    return ret;
#elif defined(GCC_X86_64)
    word ret;
    asm( "divq %4" "\n\t" "movq %%rax,%1"
    : "=d"(*r),"=g"(ret)
    :  "0"(*r), "a"(a),"g"(b)
    );
    return ret;
#else
    word qh, ql, 
         bh, bl,
         rh, rl, tp;
    ulong scale;
    if ( scale = WORDSIZE - word_log2(b) - 1 ) {
     b <<= scale;
    *r = (*r << scale) | (a >> (WORDSIZE-scale));
     a <<= scale;
    }

    bl = LO(b), bh = HI(b);
    
    qh = *r / bh;
    rh = *r - qh * bh;
    tp = qh * bl;
    rh = (rh << HWSIZE) | HI(a);
    if (rh<tp) {
        qh--,rh+=b;
        if (rh<tp && rh>=b)
            qh--,rh+=b;
    }
    rh -= tp;

    ql = rh / bh;
    rl = rh - ql * bh;
    tp = ql * bl;
    rl = (rl << HWSIZE) | LO(a);
    if (rl<tp) {
        ql--,rl+=b;
        if (rl<tp && rl>=b)
            ql--,rl+=b;
    }
    
    *r = (rl - tp) >> scale;
    return (qh<<HWSIZE)|ql;
#endif
}



/* TODO: Check if this works on *nix */
/* Calculate [b][a] / [d]. Return the quotient. */
INLINE word div_without_remainder(word b, word a, word d) {
#if defined(MVC_X86_32)
    __asm {
        mov eax,a
        mov edx,b
        div d
    }
#elif defined(GCC_X86_32)
    word ret;
    asm( "divl %3" "\n\t" "movl %%eax,%1"
    :  "=g"(ret)
    : "a"(a), "d"(b), "g"(d)
    );
    return ret;
#elif defined(GCC_X86_64)
    word ret;
    asm( "divq %3" "\n\t" "movq %%rax,%1"
    : "=g"(ret)
    : "a"(a), "d"(b), "g"(d)
    );
    return ret;
#else
    word qh, ql, 
         bh, bl,
         rh, rl, tp;
    ulong scale;
    if ( scale = WORDSIZE - word_log2(d) - 1 ) {
     d <<= scale;
     b = (b << scale) | (a >> (WORDSIZE-scale));
     a <<= scale;
    }

    bl = LO(d), bh = HI(d);
    
    qh = b / bh;
    rh = b - qh * bh;
    tp = qh * bl;
    rh = (rh << HWSIZE) | HI(a);
    if (rh<tp) {
        qh--,rh+=d;
        if (rh<tp && rh>=d)
            qh--,rh+=d;
    }
    rh -= tp;

    ql = rh / bh;
    rl = rh - ql * bh;
    tp = ql * bl;
    rl = (rl << HWSIZE) | LO(a);
    if (rl<tp) {
        ql--,rl+=d;
        if (rl<tp && rl>=d)
            ql--,rl+=d;
    }
    return (qh<<HWSIZE)|ql;
#endif
}

#else 

word div_with_remainder(word a, word b, word *r);
word div_without_remainder(word a, word b, word d);

#endif 