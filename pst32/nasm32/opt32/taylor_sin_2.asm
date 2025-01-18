section .data
    factorials dd 1.0, 6.0, 120.0, 5040.0  ; Fattoriali
    signs dd 1.0, -1.0, 1.0, -1.0          ; Segni alternati

section .text
    global taylor_sin_2

taylor_sin_2: ; x è in xmm0
    ; Calcolo dei poteri di x
    push ebp
    mov ebp, esp

    mov ecx, [ebp+8]

    movss xmm0, [ecx] ; xmm0 = [X,...,...,...]
    shufps xmm0, xmm0, 0x00 ; xmm0 = [x,x,x,x]
    
    movss xmm7, xmm0
    
    movss xmm1, xmm0           ; xmm1 = [x, x, x, x]
    mulss xmm1, xmm0            ; xmm1 = x^2 (secondo termine) xmm1 = [x^2,x^2,x^2,x^2]

    movss xmm2, xmm1           ; xmm2 = x^2
    mulss xmm2, xmm0            ; xmm2 = x^3

    movss xmm3, xmm2           ; xmm3 = x^3
    mulss xmm3, xmm1            ; xmm3 = x^5*4 volte

    movss xmm4, xmm3           ; xmm4 = x^5
    mulss xmm4, xmm1            ; xmm4 = x^7

    ; Organizza x, x^3, x^5, x^7 in xmm0
    movss xmm0, xmm2           ; xmm0 = [x^3, x^3, x^3, x^3]
    shufps xmm0, xmm3, 0x88     ; xmm0 = [x^3, x^5, x^3, x^5]
    shufps xmm0, xmm4, 0xCC     ; xmm0 = [x^3, x^5, x^7, x^7]
    shufps xmm0, xmm0, 0x1B; xmm0 = [x^7, x^3, x^5, x^7]

    ; Correggi il primo elemento per contenere x
    movss xmm0, xmm7       ; xmm0[0] = x (argomento passato alla funzione)

    ; Divisione per i fattoriali
    movaps xmm5, [factorials]   ; xmm5 = [1.0, 6.0, 120.0, 5040.0]
    divps xmm0, xmm5            ; xmm0 = [x/1!, x^3/3!, x^5/5!, x^7/7!]

    ; Applicazione dei segni
    movaps xmm6, [signs]        ; xmm6 = [1.0, -1.0, 1.0, -1.0]
    mulps xmm0, xmm6            ; xmm0 *= [1.0, -1.0, 1.0, -1.0]

    ; Somma orizzontale dei termini
    haddps xmm0, xmm0           ; xmm0 = [term0 + term1, term2 + term3, ...]
    haddps xmm0, xmm0           ; xmm0 = [result, ..., ...]

    movss [ecx], xmm0
    mov esp, ebp
    pop ebp

    ret
