section .text
    global prod_mat_opt

    prod_mat_opt:
        push ebp
        mov ebp, esp

        mov ecx, [ebp+8] ; puntatore a dist
        mov edx, [ebp+12]; puntatore a newv

        movss xmm0, [ecx]
        shufps xmm0, xmm0, 0x00; xmm0 = [d,d,d,d]

        movaps xmm1, [edx]; carico i primi 4 di rot in xmm1
        mulps xmm1, xmm0
        movaps [edx], xmm1

        mov esp, ebp
        pop ebp

        ret



