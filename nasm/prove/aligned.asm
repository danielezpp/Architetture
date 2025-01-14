section .text
    global process_array

process_array:
    push ebp
    mov ebp, esp

    push ebx
    push esi

    mov ebx, [ebp+8] ;puntatore all'array
    mov esi, [ebp+12] ;n

    .loop:
        cmp esi, 0
        je .end_loop

        movaps xmm0, [ebx] 
        addps xmm0, xmm0
        movaps [ebx], xmm0

        sub esi, 4
        add ebx, 16
        jmp .loop
    
    .end_loop:
        pop esi
        pop ebx
        
        mov esp, ebp
        pop ebp
        ret