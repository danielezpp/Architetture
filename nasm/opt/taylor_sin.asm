section .data
    factorials dd 1.0, 6.0, 120.0, 5040.0
    signs dd -1.0, 1.0, -1.0, 1.0

section .text
    global taylor_sin

    taylor_sin: ;x è in xmm0
        movaps xmm1, [factorials]
        movaps xmm2, [signs]

        movss xmm3,xmm0           ; Copia x in xmm3
        mulss xmm3, xmm0           ; xmm3 = x^2
        movss xmm4, xmm3          ; Copia x^2 in xmm4
        mulss xmm4, xmm0           ; xmm4 = x^3
        movss xmm5, xmm4          ; Copia x^3 in xmm5
        mulss xmm5, xmm3           ; xmm5 = x^5
        movss xmm6, xmm5          ; Copia x^5 in xmm6
        mulss xmm6, xmm3           ; xmm6 = x^7 

        pslld xmm0, xmm0
        movss xmm0, xmm4
        pslld xmm0, xmm0
        movss xmm0, xmm5
        pslld xmm0, xmm0
        movss xmm0, xmm6

        divps xmm0, xmm1
        mulps xmm0, xmm2

        haddps xmm0, xmm0
        haddps xmm0, xmm0

        ret



