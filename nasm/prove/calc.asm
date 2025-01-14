section .data
    two: dd 2.0          ; Definisce il valore 2.0

section .text
    global process       ; Funzione accessibile dal C

process:
    ; La funzione prende un puntatore all'array v e un intero n
    push ebp
    mov ebp, esp

    push ebx
    push esi
    push edi

    mov ecx, [ebp+8]     ; Carica il primo parametro (puntatore v) in ecx
    mov esi, [ebp+12]    ; Carica il secondo parametro (n) in esi

.loop:
    cmp esi, 0           ; Controlla se n è zero (fine array)
    je .end_loop

    ; Carica il valore di v[i] in xmm0
    movss xmm0, [ecx]    
    ; Aggiungi 2.0 al valore in xmm0
    addss xmm0, [two]    
    ; Salva il risultato di xmm0 indietro in v[i]
    movss [ecx], xmm0    

    ; Sposta al prossimo elemento dell'array
    add ecx, 4           
    dec esi               ; Decrementa n
    jmp .loop             ; Ripeti il ciclo

.end_loop:
    pop edi               ; Ripristina i registri
    pop esi
    pop ebx
    mov esp, ebp
    pop ebp
    ret                   ; Torna al chiamante
