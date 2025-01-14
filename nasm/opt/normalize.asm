section .text 
    global normalize_v

    normalize_v:
        push ebp
        mov ebp, esp

        mov ecx, [ebp+8]

        movaps xmm0, [ecx]
        movaps xmm1, xmm0
        mulps xmm1, xmm1
        haddps xmm1, xmm1
        haddps xmm1, xmm1
        sqrtss xmm1, xmm1 ; magnitudine in xmm1[0]
        shufps xmm1, xmm1, 0x00
        
        divps xmm0, xmm1
        movaps [ecx], xmm0

        mov esp, ebp
        pop ebp

        ret

