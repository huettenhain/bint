# Big Integer Library 

This was a pet project during my undergraduate studies, I maintain it here
purely out of nostalgia. The code was written almost 10 years ago and I would 
do many things differently today.

If you are looking for a good big integer library in C, use the 
[GMP](https://gmplib.org/) instead, seriously. However, if you are interested
in big integer arithmetics, feel free to browse through this barely documented
code. If I find the time, I will try to improve the code quality and maybe
insert some comments, or even documentation.

# Fast Fourier Multiplication 

Most of the algorithms are pretty standard, the exception being the FFT 
multiplication. For the x86 architecture, I implemented the FFT in pure 
assembler because it allowed for a manual optimization in the FPU computations
that beat the C compiler at that time. I have not tested it, but I do not
believe that I can beat a C compiler on x64 with the fancy new FPU 
instructions, so I implemented an FFT in C for x64. 



