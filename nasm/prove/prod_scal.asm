section .text
    global prod_scal

    prod_scal:
        push ebp
        mov ebp, esp

        push esi

        mov ecx, [ebp+8] ;v1
        mov edx, [ebp+12] ;v2
        mov esi, [ebp+16] ;n

        xorps xmm2, xmm2

        .loop:
            cmp esi, 0
            je .end_loop

            movaps xmm0, [ecx]
            movaps xmm1, [edx]
            mulps xmm0, xmm1
            haddps xmm0, xmm0
            haddps xmm0, xmm0
            addss xmm2, xmm0

            add ecx, 16
            add edx, 16
            sub esi, 4
            jmp .loop
        .end_loop:
        movss [esp-4], xmm2   ; Salva il risultato dal registro XMM in memoria temporanea
        fld dword [esp-4]     ; Carica il risultato in st(0), il registro della FPU

            pop esi
            mov esp, ebp
            pop ebp
            ret

            
        