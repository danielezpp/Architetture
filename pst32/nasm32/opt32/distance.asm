section .text
    global get_distance

    get_distance:
        push ebp
        mov ebp, esp

        mov ecx, [ebp+8]  ; puntatore a v
        mov edx, [ebp+12] ; puntatore a w

        movaps xmm6, [ecx]
        movaps xmm7, [edx]
        subps xmm6, xmm7
        mulps xmm6, xmm6
        haddps xmm6, xmm6
        haddps xmm6, xmm6
        sqrtps xmm6, xmm6

        sub esp, 4
        movss [esp], xmm6
        fld dword [esp]
        add esp, 4

        mov esp, ebp
        pop ebp
        ret

        