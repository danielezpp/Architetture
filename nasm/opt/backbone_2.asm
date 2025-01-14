section .text
    global backbone_2

    backbone_2:
        push ebp
        mov ebp, esp
        ;push esi

        mov ecx, [ebp+8] ; puntatore a v1
        mov edx, [ebp+12]; puntatore a coords
        mov eax, [ebp+16]; eax = index

        ;movaps xmm0, [ecx]; primi 4 di v1 in xmm0
        ;shufps xmm0, xmm0, 0x00; xmm0 = [d,d,d,d]

        lea esi, [edx+eax*4-4*3]
        movups xmm1, [esi]; carico coords da [index-3] a [index] in xmm1
        lea esi, [edx+eax*4+4*3]
        movups xmm2, [esi]; carico coords da [index-6] a [index -3] in xmm2

        subps xmm1, xmm2
        movups [ecx], xmm1

        ;pop esi
        mov esp, ebp
        pop ebp

        ret

