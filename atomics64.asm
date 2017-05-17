PUBLIC div_with_remainder
PUBLIC div_without_remainder

.data
.code

; x64 calling convention passes RCX, RDX, R8, R9
; and expects return value in RAX.

div_with_remainder PROC
 MOV RAX, RCX
 MOV R10, RDX
 MOV RDX, QWORD PTR [R8]
 DIV R10
 MOV QWORD PTR [R8], RDX
 RET
div_with_remainder ENDP

div_without_remainder PROC
 MOV RAX, RDX
 MOV RDX, RCX
 DIV R8
 RET
div_without_remainder ENDP

END