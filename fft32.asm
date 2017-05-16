PUBLIC _fft

.486
.model flat

.data
.code
        
_fft PROC
; ARGUMENTS:
; [-------] ESP
; [RETADDR]     [ESP+0]
; [POINTER] f = [ESP+4]
; [INTEGER] n = [ESP+8]
; [INTEGER] s = [ESP+12]
;
; no stackframe. I never understood who needs those.

		mov  ebx, [ESP+8]       ; ebx contains length of function
        mov  edx, ebx           ; so does edx

        finit                   ; initialize Floating Point Unit (not?)
        fstenv [esp-28]

        fild dword ptr [ESP+12] ; load forward/backwards indicator
        fld st(0)               ; push another copy of s on the stack
        fabs                    ; |s|,s on stack
        fdivp st(1),st(0)       ; set ST(0) to s/|s| = sgn(s)
        fist dword ptr [ESP+12] ; store the result in s, do not pop
        fld1                    ; load 1 onto stack
        fldpi                   ; load pi onto stack
        fscale                  ; scaling pi by 2
        fstp st(1)              ; now, ST(0) = 2*pi
        fmulp st(1),st(0)       ; ST(0) = 2*pi*sgn(s)
        fidiv dword ptr [ESP+8] ; ST(0) = 2*pi*sng(s)/n 
        fsincos                 ; cos(ST(0)),sin(ST(0)) on stack.

    ;  We now have w on stack:
    ;  ST(0) = Re(w)
    ;  ST(1) = Im(w)

        mov  esi, [ESP+4]            ; esi = f
        cmp  dword ptr [ESP+12], 1   ; check,
        je   START                   ; do not scale if s != 1
        fild dword ptr [ESP+8]       ; push n onto FPU stack
        mov  ecx, ebx                ; n values to be scaled  
        shl  ecx, 1                  ; each value consists of 2 doubles
SCALE:  fld  qword ptr [esi+8*ecx-8] ; push double on the stack
        fdiv ST(0),ST(1)             ; scale by n
        fstp qword ptr [esi+8*ecx-8] ; store scaled double 
        loop SCALE                   ; loop until ecx=0
       fcomp                         ; pop n off the FPU stack again.

    ;  In terms of the C version of the algorithm, we will use
    ;  the ebx register as N, the ecx register as j, the eax
    ;  register as k and the edx register as M.
 
START:  cmp  ebx,1                   ; Only contiune to transform as 
        jng  ENDFUN                  ; long as N > 1.
        shr  edx,1                   ; Set M := N >> 1;
        xor  ecx,ecx                 ; start new loop with j=0
INNER1: fldz                         ; push 0 onto FPU stack
        fld1                         ; push 1 onto FPU stack

    ;  we now have w and the initial z on the FPU stack:
    ;  ST(0) = Re(z)
    ;  ST(1) = Im(z)
    ;  ST(2) = Re(w)
    ;  ST(3) = Im(w)

        xor  eax,eax                 ; start new loop with k=0
INNER2: add  eax,ecx                 ; for indexing, we need k+j
        mov  edi,eax                 ; and
        add  edi,edx                 ; k+j+M

        shl  eax,4                   ; we are indexing 16-byte COMPLEX
        shl  edi,4                   ; values, scale indices

        fld  qword ptr [esi+edi+8]   ; ST(2)  = f[j+k+M].imag
        fld  qword ptr [esi+eax+8]   ; ST(1)  = f[j+k].imag
        fld  ST(1)                   ; ST(0)  = f[j+k+M].imag
        fadd ST(0),ST(1)             ; ST(0)  = f[j+k].imag + f[j+k+M].imag
        fstp qword ptr [esi+eax+8]   ; f[j+k] = f[j+k].imag + f[j+k+M].imag
        fsubrp ST(1),ST(0)           ; ST(0)  = f[j+k].imag - f[j+k+M].imag

    ;  Now, the same procedure is performed for the real parts of
    ;  the two complex variables:

        fld  qword ptr [esi+edi]
        fld  qword ptr [esi+eax]
        fld  ST(1)
        fadd ST(0),ST(1)
        fstp qword ptr [esi+eax]
        fsubrp ST(1),ST(0)

    ;  the FPU stack now looks like this:
    ;
    ;  ST(0) = real part of f[j+k]-f[j+k+M] =: a
    ;  ST(1) = imag part of a
    ;  ST(2) = real of z
    ;  ST(3) = imag of z
    ;  ST(4) = real of w
    ;  ST(5) = imag of w
    ;
    ;  we want to calculate a*z and store it in f[j+k+M].

        fld   ST(0)                  ; push Re(a)
        fmul  ST(0),ST(3)            ; ST(0) <- Re(z)*Re(a)
        fld   ST(2)                  ; push Im(a)
        fmul  ST(0),ST(5)            ; ST(0) <- Im(a)*Im(z)
        fsubp ST(1),ST(0)            ; TOS   <- ST(1)-ST(0)
        fstp  qword ptr [esi+edi]    ; pop to f[j+k+M].real
    ; ----------------------------- ; stack back to [a,z,w]
        fmul  ST(0),ST(3)            ; ST(0) <- Re(a)*Im(z)
        fxch  ST(1)                  ; ST(0) <-> ST(1)
        fmul  ST(0),ST(2)            ; ST(1) <- Im(a)*Re(z)
        faddp ST(1),ST(0)            ; TOS   <- ST(0)+ST(1)
        fstp  qword ptr [esi+edi+8]  ; pop to f[j+k+M].imag
    
    ;  We now have only z and w on the FPU stack and for the 
    ;  next iteration, we want to set z := z*w.

        fld   ST(0)                  ; push Re(z)
        fmul  ST(0),ST(4)            ; ST(0) <- Re(z)*Im(w)
        fld   ST(2)                  ; push Im(z)
        fmul  ST(0),ST(4)            ; ST(0) <- Im(z)*Re(w)
        faddp ST(1),ST(0)            ; TOS   <- Im(z*w)
        fxch  ST(2)                  ; set new Im(z)
        fmul  ST(0),ST(4)            ; ST(0) <- Im(z)*Im(w)
        fxch  ST(1)                  ; ST(0) <-> ST(1)
        fmul  ST(0),ST(3)            ; ST(0) <- Re(z)*Re(w)
        fsubp ST(1),ST(0)            ; ST(0) <- -Re(z*w) 
        fchs                         ; ST(0) <- Re(z*w)

        shr   eax,4                  ; unscale indices again

        sub   eax,ecx                ; eax := (k+j)-j = k
        inc   eax                    ; increase by one
        cmp   eax,edx                ; check whether k<M
        jl    INNER2                 ; if yes, loop again
                                     ; the inner loop ends here:
        fcompp                       ; pop z off the stack
        add   ecx,ebx                ; j += N
        cmp   ecx,[ESP+8]            ; check whether j < n
        jl    INNER1                 ; if yes, loop again

    ;  the outer loop ends here. We now want to calculate 
    ;  w := w*w and simply jump back to START where the 
    ;  N>1 condition is checked.

        fld   ST(1)                  ; Push Im(w) onto stack
        fmul  ST(0),ST(1)            ; ST(0) := Re(w)*Im(w)
        fld   ST(0)                  ; Push ST(0) again
        faddp ST(1),ST(0)            ; ST(0) = 2*Re(w)*Im(w)
        fxch  ST(2)                  ; exchange it with Im(w)
        fmul  ST(0),ST(0)            ; ST(0) := Im(w)*Im(w)
        fxch  ST(1)                  ; exchange it with Re(w)
        fmul  ST(0),ST(0)            ; ST(0) := Re(w)*Re(w)
        fsubp ST(1),ST(0)            ; ST(0) = Im(w)*Im(w)-Re(w)*Re(w)
        fchs                         ; invert sign, done

        mov   ebx,edx                ; N := M
        jmp   START                  ; and again

    ;  The transformation has now been performed completely. All
    ;  that is left to be done is doing the bit reversal.

ENDFUN: xor   edx,edx                ; edx will be first index
        mov   eax,edx                ; eax the second
        mov   ebx, dword ptr [ESP+8] ; ebx stores n
BREV:   cmp   edx,eax                ; only swap values if 
        jle   NOSWAP                 ; edx > eax
        shl   edx,4                  
        shl   eax,4
        fld   qword ptr [esi+edx]
        fld   qword ptr [esi+edx+8]
        fld   qword ptr [esi+eax]
        fld   qword ptr [esi+eax+8]
        fstp  qword ptr [esi+edx+8]
        fstp  qword ptr [esi+edx]
        fstp  qword ptr [esi+eax+8]
        fstp  qword ptr [esi+eax]
        shr   edx,4
        shr   eax,4
NOSWAP: mov   ecx,ebx                ; calculate bitreversal:
CLOOP:  shr   ecx,1                  ; ecx is of the form 1<<k.
        xor   edx,ecx                ; XOR ecx to edx until the result of
        cmp   edx,ecx                ; that XOR is not 0, hence ecx<=edx
        jl    CLOOP                  ; otherwise we have a 'carry'.
        inc   eax                    ; increase normal index
        cmp   eax,ebx                ; and loop until we're done
        jl    BREV
        mov   eax,esi                ; return f;

        fldenv [esp-28]

		ret

_fft ENDP


END