dnl
dnl Copyright (C) 2023, 2024 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

define(`rp', `%rdi')
define(`ap', `%rsi')

define(`s0', `%rcx')
define(`s1', `%r8')
define(`s2', `%r9')
define(`s3', `%r10')
define(`s4', `%r11')
define(`s5', `%rbx')
define(`s6', `%rbp')
define(`s7', `%r12')
define(`s8', `%r13')
define(`s9', `%r14')
define(`s10', `%r15')

define(`sx', `%rax')

	TEXT

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_1)
	mov	0*8(ap), %rdx
	mulx	%rdx, s0, sx
	mov	s0, 0*8(rp)
	mov	sx, 1*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_2)
	mov	0*8(ap), %rdx
	mulx	1*8(ap), s1, s2		C a0 a1
	mulx	%rdx, s3, s4		C a0^2
	add	s1, s1
	adc	s2, s2
	mov	1*8(ap), %rdx
	mulx	%rdx, s0, sx		C a1^2
	mov	s3, 0*8(rp)
	adc	$0, sx
	add	s4, s1
	adc	s0, s2
	adc	$0, sx
	mov	s1, 1*8(rp)
	mov	s2, 2*8(rp)
	mov	sx, 3*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_3)
	mov	0*8(ap), %rdx
	mulx	1*8(ap), sx, s0		C a0 a1
	mulx	2*8(ap), s1, s2		C a0 a2

	add	s0, s1
	mov	1*8(ap), %rdx
	mulx	2*8(ap), s3, s4		C a1 a2

	mov	$0, R32(s0)
	adc	s2, s3
	adc	s0, s4

	add	sx, sx
	adc	s1, s1
	adc	s3, s3
	adc	s4, s4

	mov	0*8(ap), %rdx
	mulx	%rdx, %rdx, s2		C a0^2
	adc	R32(s0), R32(s0)
	mov	%rdx, 0*8(rp)
	add	s2, sx
	mov	sx, 1*8(rp)

	mov	1*8(ap), %rdx
	mulx	%rdx, s2, sx		C a1^2
	adc	s2, s1
	adc	sx, s3

	mov	2*8(ap), %rdx
	mulx	%rdx, s2, sx		C a2^2
	mov	s1, 2*8(rp)
	mov	s3, 3*8(rp)
	adc	s4, s2
	adc	s0, sx
	mov	s2, 4*8(rp)
	mov	sx, 5*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_4)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	xor	R32(s0), R32(s0)

	mulx	1*8(ap), s1, s4		C a0 a1
	mulx	2*8(ap), s3, s2		C a0 a2
	mulx	3*8(ap), s5, s6		C a0 a3
	mov	1*8(ap), %rdx
	adox	s4, s3
	adox	s2, s5
	mulx	2*8(ap), s4, s2		C a1 a2
	adcx	s4, s5
	adcx	s2, s6
	mulx	3*8(ap), s4, s2		C a1 a3
	mov	2*8(ap), %rdx
	adox	s4, s6
	adox	s0, s2
	mulx	3*8(ap), %rdx, s4	C a1 a3
	adc	%rdx, s2
	adc	s0, s4
	mov	0*8(ap), %rdx
	add	s1, s1
	adc	s3, s3
	adc	s5, s5
	adc	s6, s6
	adc	s2, s2
	adc	s4, s4
	setc	R8(s0)

	mulx	%rdx, %rdx, sx		C a0^2
	mov	%rdx, 0*8(rp)
	add	sx, s1
	mov	1*8(ap), %rdx
	mov	s1, 1*8(rp)
	mulx	%rdx, %rdx, sx		C a1^2
	adc	%rdx, s3
	adc	sx, s5
	mov	2*8(ap), %rdx
	mov	s3, 2*8(rp)
	mov	s5, 3*8(rp)
	mulx	%rdx, s1, sx		C a2^2
	adc	s1, s6
	adc	sx, s2
	mov	3*8(ap), %rdx
	mov	s6, 4*8(rp)
	mov	s2, 5*8(rp)
	mulx	%rdx, s1, sx		C a3^2
	adc	s1, s4
	adc	s0, sx
	mov	s4, 6*8(rp)
	mov	sx, 7*8(rp)

	pop	s6
	pop	s5

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_5)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	push	s7
	xor	R32(s0), R32(s0)

	mulx	1*8(ap), sx, s2		C a0 a1
	mulx	2*8(ap), s3, s4		C a0 a2
	mulx	3*8(ap), s5, s1		C a0 a3
	mulx	4*8(ap), s6, s7		C a0 a4
	adcx	s2, s3
	adcx	s4, s5
	adcx	s1, s6
	C x, 3, 5, 6, 7

	mov	1*8(ap), %rdx
	mulx	2*8(ap), s2, s4		C a1 a2
	adcx	s0, s7
	adcx	s2, s5
	mulx	3*8(ap), s2, s1		C a1 a3
	adcx	s4, s6
	adox	s2, s6
	mulx	4*8(ap), s2, s4		C a1 a4
	adox	s1, s7
	adcx	s2, s7
	adcx	s0, s4
	C x, 3, 5, 6, 7, 4

	mov	2*8(ap), %rdx
	mulx	3*8(ap), s1, s2		C a2 a3
	adcx	s1, s7
	adcx	s2, s4
	mulx	4*8(ap), s1, s2		C a2 a4
	adox	s1, s4
	C x, 3, 5, 6, 7, 4, 2

	mov	3*8(ap), %rdx
	mulx	4*8(ap), %rdx, s1	C a3 a4
	adox	s0, s2
	adc	%rdx, s2
	adc	s0, s1
	C x, 3, 5, 6, 7, 4, 2, 1

	mov	0*8(ap), %rdx
	add	sx, sx
	adc	s3, s3
	mov	sx, -1*8(%rsp)
	adc	s5, s5
	adc	s6, s6
	mulx	%rdx, %rdx, sx		C a0^2
	adc	s7, s7
	adc	s4, s4
	adc	s2, s2
	adc	s1, s1
	setc	R8(s0)
	C x, 3, 5, 6, 7, 4, 2, 1, 0

	mov	%rdx, 0*8(rp)
	add	-1*8(%rsp), sx
	mov	sx, 1*8(rp)

	mov	1*8(ap), %rdx
	mulx	%rdx, %rdx, sx		C a1^2
	adc	%rdx, s3
	adc	sx, s5
	mov	s3, 2*8(rp)

	mov	2*8(ap), %rdx
	mulx	%rdx, s3, sx		C a2^2
	mov	s5, 3*8(rp)
	adc	s3, s6
	adc	sx, s7
	mov	s6, 4*8(rp)

	mov	3*8(ap), %rdx
	mulx	%rdx, s3, sx		C a3^2
	mov	s7, 5*8(rp)
	adc	s3, s4
	adc	sx, s2
	mov	s4, 6*8(rp)

	mov	4*8(ap), %rdx
	mulx	%rdx, s3, sx		C a4^2
	mov	s2, 7*8(rp)
	adc	s1, s3
	adc	s0, sx
	mov	s3, 8*8(rp)
	mov	sx, 9*8(rp)

	pop	s7
	pop	s6
	pop	s5

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_6)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	push	s7
	push	s8
	push	s9
	xor	R32(s0), R32(s0)

	mulx	1*8(ap), sx, s2		C a0 a1
	mulx	2*8(ap), s3, s4		C a0 a2
	mulx	3*8(ap), s5, s1		C a0 a3
	mulx	4*8(ap), s6, s9		C a0 a4
	mulx	5*8(ap), s7, s8		C a0 a5
	adox	s2, s3
	adox	s4, s5
	adox	s1, s6
	adox	s9, s7
	adox	s0, s8
	C x, 3, 5, 6, 7, 8

	mov	1*8(ap), %rdx
	mulx	2*8(ap), s2, s4		C a1 a2
	mulx	3*8(ap), s1, s9		C a1 a3
	adcx	s2, s5
	adcx	s4, s6
	adox	s1, s6
	adox	s9, s7
	mulx	4*8(ap), s2, s4		C a1 a4
	mulx	5*8(ap), s1, s9		C a1 a5
	adcx	s2, s7
	adcx	s4, s8
	adox	s1, s8
	C x, 3, 5, 6, 7, 8, 9

	mov	2*8(ap), %rdx
	mulx	3*8(ap), s1, s2		C a2 a3
	adcx	s0, s9
	adox	s0, s9
	adcx	s1, s7
	adcx	s2, s8
	mulx	4*8(ap), s4, s1		C a2 a4
	mulx	5*8(ap), %rdx, s2	C a2 a5
	adox	s4, s8
	adox	s1, s9
	adcx	%rdx, s9
	C x, 3, 5, 6, 7, 8, 9, 2

	mov	3*8(ap), %rdx
	mulx	4*8(ap), s1, s4		C a3 a4
	adcx	s0, s2
	adcx	s1, s9
	mulx	5*8(ap), %rdx, s1	C a3 a5
	adcx	s4, s2
	adox	%rdx, s2
	adox	s0, s1
	C x, 3, 5, 6, 7, 8, 9, 2, 1

	mov	4*8(ap), %rdx
	mulx	5*8(ap), %rdx, s4	C a4 a5
	adcx	%rdx, s1
	adcx	s0, s4
	C x, 3, 5, 6, 7, 8, 9, 2, 1, 4

	add	sx, sx
	adc	s3, s3
	mov	0*8(ap), %rdx
	mov	sx, -1*8(%rsp)
	adc	s5, s5
	adc	s6, s6
	mulx	%rdx, %rdx, sx
	adc	s7, s7
	adc	s8, s8
	adc	s9, s9
	adc	s2, s2
	adc	s1, s1
	adc	s4, s4
	setc	R8(s0)
	C x, 3, 5, 6, 7, 8, 9, 2, 1, 4, 0

	mov	%rdx, 0*8(rp)
	add	-1*8(%rsp), sx
	mov	sx, 1*8(rp)

	mov	1*8(ap), %rdx
	mulx	%rdx, %rdx, sx		C a1^2
	adc	%rdx, s3
	adc	sx, s5
	mov	s3, 2*8(rp)

	mov	2*8(ap), %rdx
	mulx	%rdx, s3, sx		C a2^2
	mov	s5, 3*8(rp)
	adc	s3, s6
	adc	sx, s7
	mov	s6, 4*8(rp)

	mov	3*8(ap), %rdx
	mulx	%rdx, s3, sx		C a3^2
	mov	s7, 5*8(rp)
	adc	s3, s8
	adc	sx, s9
	mov	s8, 6*8(rp)

	mov	4*8(ap), %rdx
	mulx	%rdx, s3, sx		C a4^2
	mov	s9, 7*8(rp)
	adc	s3, s2
	adc	sx, s1
	mov	s2, 8*8(rp)

	mov	5*8(ap), %rdx
	mulx	%rdx, s3, sx		C a5^2
	mov	s1, 9*8(rp)
	adc	s4, s3
	adc	s0, sx
	mov	s3, 10*8(rp)
	mov	sx, 11*8(rp)

	pop	s9
	pop	s8
	pop	s7
	pop	s6
	pop	s5

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_7)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	push	s7
	push	s8
	push	s9
	push	s10
	xor	R32(s0), R32(s0)

	mulx	1*8(ap), sx, s1		C a0 a1
	mulx	2*8(ap), s3, s2		C a0 a2
	mulx	3*8(ap), s5, s4		C a0 a3
	mulx	4*8(ap), s6, s10	C a0 a4
	mov	sx, -2*8(%rsp)
	adox	s1, s3
	adox	s2, s5
	mov	s3, -1*8(%rsp)
	mulx	5*8(ap), s7, s1		C a0 a5
	mulx	6*8(ap), s8, s9		C a0 a6
	adox	s4, s6
	adox	s10, s7
	adox	s1, s8
	adox	s0, s9
	C (-, -,) 5, 6, 7, 8, 9
	C x, 1, 2, 3, 4, 10

	mov	1*8(ap), %rdx
	mulx	2*8(ap), sx, s1		C a1 a2
	mulx	3*8(ap), s2, s3		C a1 a3
	adox	sx, s5
	adox	s1, s6
	adcx	s2, s6
	adcx	s3, s7
	mulx	4*8(ap), sx, s1		C a1 a4
	mulx	5*8(ap), s2, s3		C a1 a5
	mulx	6*8(ap), s4, s10	C a1 a6
	adox	sx, s7
	adox	s1, s8
	adcx	s2, s8
	adcx	s3, s9
	adox	s4, s9
	adox	s0, s10
	adcx	s0, s10
	C (-, -,) 5, 6, 7, 8, 9, 10
	C x, 1, 2, 3, 4

	mov	2*8(ap), %rdx
	mulx	3*8(ap), sx, s1		C a2 a3
	mulx	4*8(ap), s2, s3		C a2 a4
	adox	sx, s7
	adox	s1, s8
	adcx	s2, s8
	adcx	s3, s9
	mulx	5*8(ap), sx, s1		C a2 a5
	mulx	6*8(ap), s2, s3		C a2 a6
	adox	sx, s9
	adox	s1, s10
	adcx	s2, s10
	adcx	s0, s3
	adox	s0, s3
	C (-, -,) 5, 6, 7, 8, 9, 10, 3
	C x, 1, 2, 4

	mov	3*8(ap), %rdx
	mulx	4*8(ap), sx, s1		C a3 a4
	mulx	5*8(ap), s2, s4		C a3 a5
	adox	sx, s9
	adox	s1, s10
	mulx	6*8(ap), sx, s1		C a3 a6
	adcx	s2, s10
	adcx	s4, s3
	adox	sx, s3
	adox	s0, s1
	C (-, -,) 5, 6, 7, 8, 9, 10, 3, 1
	C x, 2, 4

	mov	4*8(ap), %rdx
	mulx	5*8(ap), sx, s2		C a4 a5
	mulx	6*8(ap), %rdx, s4	C a4 a6
	adox	sx, s3
	adox	s2, s1
	adcx	%rdx, s1
	adox	s0, s4
	C (-, -,) 5, 6, 7, 8, 9, 10, 3, 1, 4
	C x, 2

	mov	5*8(ap), %rdx
	mulx	6*8(ap), sx, s2		C a5 a6
	adc	sx, s4
	adc	s0, s2
	C (-, -,) 5, 6, 7, 8, 9, 10, 3, 1, 4, 2
	C x

	shlq	-2*8(%rsp)
	mov	-1*8(%rsp), sx
	mov	0*8(ap), %rdx
	adc	sx, sx
	adc	s5, s5
	mov	sx, -1*8(%rsp)
	adc	s6, s6
	adc	s7, s7
	mulx	%rdx, %rdx, sx
	adc	s8, s8
	adc	s9, s9
	adc	s10, s10
	adc	s3, s3
	adc	s1, s1
	adc	s4, s4
	adc	s2, s2
	setc	R8(s0)
	C (-, -,) 5, 6, 7, 8, 9, 10, 3, 1, 4, 2, 0
	C x

	add	-2*8(%rsp), sx
	mov	%rdx, 0*8(rp)
	mov	sx, 1*8(rp)

	mov	1*8(ap), %rdx
	mulx	%rdx, %rdx, sx		C a1^2
	adc	-1*8(%rsp), %rdx
	mov	%rdx, 2*8(rp)
	adc	sx, s5

	mov	2*8(ap), %rdx
	mulx	%rdx, %rdx, sx		C a2^2
	mov	s5, 3*8(rp)
	adc	%rdx, s6
	adc	sx, s7
	mov	s6, 4*8(rp)

	mov	3*8(ap), %rdx
	mulx	%rdx, s5, sx		C a3^2
	mov	s7, 5*8(rp)
	adc	s5, s8
	adc	sx, s9
	mov	s8, 6*8(rp)

	mov	4*8(ap), %rdx
	mulx	%rdx, s5, sx		C a4^2
	mov	s9, 7*8(rp)
	adc	s5, s10
	adc	sx, s3
	mov	s10, 8*8(rp)

	mov	5*8(ap), %rdx
	mulx	%rdx, s5, sx		C a5^2
	mov	s3, 9*8(rp)
	adc	s5, s1
	adc	sx, s4
	mov	s1, 10*8(rp)

	mov	6*8(ap), %rdx
	mulx	%rdx, s3, sx		C a6^2
	mov	s4, 11*8(rp)
	adc	s2, s3
	adc	s0, sx
	mov	s3, 12*8(rp)
	mov	sx, 13*8(rp)

	pop	s10
	pop	s9
	pop	s8
	pop	s7
	pop	s6
	pop	s5

	ret
EPILOGUE()

dnl     <<  >>
dnl -1      -x	<- Must be 0
dnl  0   x-/-x
dnl  1   x-/-x
dnl  2   x-/-x
dnl  3   x-/-x
dnl  4   x-/-x
dnl  5   x-/-x
dnl  6   x-/-x
dnl  7   x-/-x
dnl  8   x-/-x
dnl  9   x-/-x
dnl 10   x-/-x
dnl 11   x-/	<- s2
dnl 12		<- s3
dnl 13		<- s1

define(`m',`eval($1-12)*8(%rsp)')

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_8)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	push	s7

	mulx	1*8(ap), sx, s0		C a0 a1
	mulx	2*8(ap), s1, s2		C a0 a2
	mulx	3*8(ap), s3, s4		C a0 a3
	mov	sx, m(0)
	add	s0, s1
	adc	s2, s3
	mov	s1, m(1)
	mulx	4*8(ap), s5, s6		C a0 a4
	mulx	5*8(ap), sx, s0		C a0 a5
	mulx	6*8(ap), s2, s1		C a0 a6
	adc	s4, s5
	adc	s6, sx
	mulx	7*8(ap), s4, s7		C a0 a7
	adc	s0, s2
	adc	s1, s4
	adc	$0, s7
	C m(0:1), 3, 5, x, 2, 4, 7
	C 1, 6

	xor	R32(s0), R32(s0)
	mov	1*8(ap), %rdx
	mulx	2*8(ap), s1, s6		C a1 a2
	mov	s0, m(-1)
	adox	s1, s3
	adox	s6, s5
	mulx	3*8(ap), s1, s6		C a1 a3
	mov	s3, m(2)
	adcx	s1, s5
	adcx	s6, sx
	mov	s5, m(3)
	mulx	4*8(ap), s1, s6		C a1 a4
	mulx	5*8(ap), s3, s5		C a1 a5
	adox	s1, sx
	adox	s6, s2
	adcx	s3, s2
	adcx	s5, s4
	mulx	6*8(ap), s1, s6		C a1 a6
	mulx	7*8(ap), s3, s5		C a1 a7
	adox	s1, s4
	adox	s6, s7
	adcx	s3, s7
	adcx	s0, s5
	adox	s0, s5
	C m(0:3), x, 2, 4, 7, 5
	C 1, 3, 6

	mov	2*8(ap), %rdx
	mulx	3*8(ap), s1, s3		C a2 a3
	adcx	s1, sx
	adcx	s3, s2
	mov	sx, m(4)
	mulx	4*8(ap), s6, s1		C a2 a4
	mulx	5*8(ap), s3, sx		C a2 a5
	adox	s6, s2
	adox	s1, s4
	adcx	s3, s4
	adcx	sx, s7
	mulx	6*8(ap), s6, s1		C a2 a6
	mulx	7*8(ap), s3, sx		C a2 a7
	mov	s2, m(5)
	adox	s6, s7
	adox	s1, s5
	adcx	s3, s5
	adcx	s0, sx
	adox	s0, sx
	C m(0:5), 4, 7, 5, x
	C 1, 2, 3, 6

	mov	3*8(ap), %rdx
	mulx	4*8(ap), s1, s2		C a3 a4
	mulx	5*8(ap), s3, s6		C a3 a5
	adox	s1, s4
	adox	s2, s7
	adcx	s3, s7
	adcx	s6, s5
	mov	s4, m(6)
	mov	s7, m(7)
	mulx	6*8(ap), s1, s2		C a3 a6
	mulx	7*8(ap), s3, s4		C a3 a7
	adox	s1, s5
	adox	s2, sx
	adcx	s3, sx
	adcx	s0, s4
	adox	s0, s4
	C m(0:7), 5, x, 4
	C 1, 2, 3, 6, 7

	mov	4*8(ap), %rdx
	mulx	5*8(ap), s1, s2		C a4 a5
	mulx	6*8(ap), s3, s6		C a4 a6
	adox	s1, s5
	adox	s2, sx
	adcx	s3, sx
	adcx	s6, s4
	mulx	7*8(ap), s1, s2		C a4 a7
	mov	s5, m(8)
	mov	sx, m(9)
	adox	s1, s4
	adox	s0, s2
	C m(0:9), 4, 2
	C x, 1, 3, 5, 6, 7

	mov	5*8(ap), %rdx
	mulx	6*8(ap), sx, s1		C a5 a6
	mulx	7*8(ap), s7, s3		C a5 a7
	adox	sx, s4
	adox	s1, s2
	adcx	s7, s2
	adox	s0, s3
	mov	s4, m(10)
	mov	s2, m(11)
	C m(0:11), (2,) 3
	C x, 1, 4, 5, 6, 7

	mov	6*8(ap), %rdx
	mulx	7*8(ap), sx, s1		C a6 a7
	adc	sx, s3
	adc	s0, s1
	C m(0:11), (2,) 3, 1
	C x, 4, 5, 6, 7

	vmovdqu	m(-1), %ymm0
	vmovdqu	m(0), %ymm1
	vmovdqu	m(3), %ymm2
	vmovdqu	m(4), %ymm3
	vmovdqu	m(7), %ymm4
	vmovdqu	m(8), %ymm5

	vpsrlq	$63, %ymm0, %ymm0
	vpsllq	$1, %ymm1, %ymm1
	vpsrlq	$63, %ymm2, %ymm2
	vpsllq	$1, %ymm3, %ymm3
	vpsrlq	$63, %ymm4, %ymm4
	vpsllq	$1, %ymm5, %ymm5

	vpor	%ymm0, %ymm1, %ymm0
	vpor	%ymm2, %ymm3, %ymm2
	vpor	%ymm4, %ymm5, %ymm4

	vmovdqu	%ymm0, m(0)
	vmovdqu	%ymm2, m(4)
	vmovdqu	%ymm4, m(8)

	mov	$63, R32(s5)
	mov	$1, R32(s6)

	C s2 = sx or s2 = (s3 << 1) or (s2 >> 63)
	shrx	s5, s2, s2
	shlx	s6, s3, sx

	C s3 = s4 or s3 = (s1 << 1) or (s3 >> 63)
	shrx	s5, s3, s3
	shlx	s6, s1, s4

	C s1 =                          s1 >> 63
	shrx	s5, s1, s1

	mov	0*8(ap), %rdx
	mulx	%rdx, s5, s6		C a0^2
	or	sx, s2
	or	s4, s3
	add	m(0), s6

	mov	1*8(ap), %rdx
	mulx	%rdx, s4, sx		C a1^2
	mov	s5, 0*8(rp)
	mov	s6, 1*8(rp)
	adc	m(1), s4
	adc	m(2), sx
	mov	s4, 2*8(rp)

	mov	2*8(ap), %rdx
	mulx	%rdx, s5, s6		C a2^2
	mov	sx, 3*8(rp)
	adc	m(3), s5
	adc	m(4), s6
	mov	s5, 4*8(rp)

	mov	3*8(ap), %rdx
	mulx	%rdx, s4, sx		C a3^2
	mov	s6, 5*8(rp)
	adc	m(5), s4
	adc	m(6), sx
	mov	s4, 6*8(rp)

	mov	4*8(ap), %rdx
	mulx	%rdx, s5, s6		C a4^2
	mov	sx, 7*8(rp)
	adc	m(7), s5
	adc	m(8), s6
	mov	s5, 8*8(rp)

	mov	5*8(ap), %rdx
	mulx	%rdx, s4, sx		C a5^2
	mov	s6, 9*8(rp)
	adc	m(9), s4
	adc	m(10), sx
	mov	s4, 10*8(rp)

	mov	6*8(ap), %rdx
	mulx	%rdx, s5, s6		C a6^2
	mov	sx, 11*8(rp)
	adc	m(11), s5
	adc	s2, s6
	mov	s5, 12*8(rp)

	mov	7*8(ap), %rdx
	mulx	%rdx, s4, sx		C a7^2
	mov	s6, 13*8(rp)
	adc	s3, s4
	adc	s1, sx
	mov	s4, 14*8(rp)
	mov	sx, 15*8(rp)

	vzeroupper

	pop	s7
	pop	s6
	pop	s5

	ret
EPILOGUE()
