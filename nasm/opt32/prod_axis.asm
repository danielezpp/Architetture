section .text
    global prod_axis

    prod_axis:
        push ebp
        mov ebp, esp

        mov ecx, [ebp+8] ; puntatore ad axis
        mov edx, [ebp+12] ; puntatore a norm

        movaps xmm0, [ecx]
        movaps xmm1, xmm0
        mulps xmm0, xmm1
        haddps xmm0, xmm0
        haddps xmm0, xmm0
        divps xmm1, xmm0

        movaps [edx], xmm1; mettiamo in normalized il risultato corretto
        
        mov esp, ebp
        pop ebp

        ret