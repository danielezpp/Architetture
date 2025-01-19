section .data
    align 16
    alpha_phi dd -57.8, -57.8, -57.8, -57.8 
    alpha_psi dd -47.0, -47.0, -47.0, -47.0
    beta_phi  dd -119.0, -119.0, -119.0, -119.0
    beta_psi  dd 113.0, 113.0, 113.0, 113.0
    half      dd 0.5, 0.5, 0.5, 0.5

section .text
    global rama_energy
rama_energy:
    ; Prologo della funzione
    push ebp
    mov ebp, esp
    and esp, 0xFFFFFFF0
    push esi
    push ebx      

    xorps xmm0, xmm0             ; energy = 0.o
    mov ecx, [ebp+8]        ; ecx punta a phi
    mov edx, [ebp+12]       ; edx punta a psi
    mov ebx, [ebp+16]       ; ebx punta a energy
    mov esi, 63             ; esi è il num di iterazioni

    .loop:
        cmp esi, 0
        je .end_loop
        movaps xmm1, [ecx]; xmm1 = 4 elementi di phi
        movaps xmm3, xmm1 ;copiamo xmm1 in xmm3
        movaps xmm2, [edx]
        movaps xmm4, xmm2 ; copiamo xmm2 in xmm4
        
        subps xmm1, [alpha_phi] ; xmm1= phi[i]-alpha_phi
        mulps xmm1, xmm1    ; xmm1^2
        subps xmm2, [alpha_psi]
        mulps xmm2, xmm2
        addps xmm1, xmm2
        sqrtps xmm1, xmm1; xmm1 = alpha_dist

        subps xmm3, [beta_phi] ; xmm3= phi[i]-beta_phi
        mulps xmm3, xmm3    ; xmm3^2
        subps xmm4, [beta_psi]
        mulps xmm4, xmm4
        addps xmm3, xmm4
        sqrtps xmm3, xmm3; xmm3 = beta_dist

        minps xmm1, xmm3 ; xmm1=min(alpha_dist, beta_dist)
        mulps xmm1, [half]
        addps xmm0, xmm1 ; aggiorno somme parziali in xmm0

        dec esi
        add ecx, 16
        add edx, 16
        jmp .loop
    .end_loop:
        haddps xmm0, xmm0
        haddps xmm0, xmm0
        movss [ebx], xmm0

    pop ebx
    pop esi
    ; Epilogo della funzione
    mov esp, ebp
    pop ebp
    ret
