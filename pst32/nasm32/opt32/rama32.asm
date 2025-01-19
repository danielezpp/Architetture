section .data
    align 16
    alpha_phi dd -57.8, -57.8, -57.8, -57.8 
    alpha_psi dd -47.0, -47.0, -47.0, -47.0
    beta_phi  dd -119.0, -119.0, -119.0, -119.0
    beta_psi  dd 113.0, 113.0, 113.0, 113.0
    half      dd 0.5, 0.5, 0.5, 0.5

section .text 
    global rama_energy32

    rama_energy32:
        push ebp
        mov ebp, esp

        push edi
        push esi

        mov ecx, [ebp+8] ; ecx punta a phi
        mov edx, [ebp+12] ;edx punta a psi
        mov esi, [ebp+16] ;in esi metto q
        mov edi, [ebp+20] ;in edi metto r

        xorps xmm0, xmm0; azzero il registro che accumula energy
        xorps xmm5, xmm5; qui accumuliamo energia del ciclo r
        dec esi

        loop_q:
            cmp esi, 0
            je end_q

            movaps xmm1, [ecx]; xmm1 = 4 elementi di phi
            movaps xmm3, xmm1 ;copiamo xmm1 in xmm3
            movaps xmm2, [edx]; xmm2 = 4 elementi di psi
            movaps xmm4, xmm2 ; copiamo xmm2 in xmm4

            ;calcolo alpha_dist
            subps xmm1, [alpha_phi] ; xmm1= phi[i]-alpha_phi
            mulps xmm1, xmm1    ; xmm1^2
            subps xmm2, [alpha_psi]
            mulps xmm2, xmm2
            addps xmm1, xmm2
            sqrtps xmm1, xmm1; xmm1 = alpha_dist

            ;calcolo beta_dist
            subps xmm3, [beta_phi] ; xmm3= phi[i]-beta_phi
            mulps xmm3, xmm3    ; xmm3^2
            subps xmm4, [beta_psi]
            mulps xmm4, xmm4
            addps xmm3, xmm4
            sqrtps xmm3, xmm3; xmm3 = beta_dist

            minps xmm1, xmm3; xmm1 = min (alpha_dist, beta_dist)
            mulps xmm1, [half]
            addps xmm0, xmm1; aggiornamento somme parziali

            dec esi
            add ecx, 16
            add edx, 16
            jmp loop_q
        end_q:
            haddps xmm0, xmm0
            haddps xmm0, xmm0; in xmm0 abbiamo rama energy del ciclo quoziente
        loop_r:
            cmp edi, 0
            je end_r

            movss xmm6, [ecx]; sposto phi[i] in xmm6 e xmm7
            movss xmm7, xmm6
            movss xmm1, [edx]; xmm1 = phi[i] e xmm2 = phi[i]
            movss xmm2, xmm1
            
            ;calcolo alpha dist
            subss xmm6, [alpha_phi]
            mulss xmm6, xmm6
            subss xmm1, [alpha_psi]
            addss xmm6, xmm1
            sqrtss xmm6, xmm6 ; xmm6 = alpha_dist

            ;calcolo beta_dist
            subss xmm7, [beta_phi]
            mulss xmm7, xmm7 
            subss xmm2, [beta_psi]
            addss xmm7, xmm2
            sqrtss xmm7, xmm2 ; xmm7 = beta_dist

            minss xmm6, xmm7; min tra alpha e beta dist
            mulss xmm6, [half]
            addss xmm5, xmm6

            dec edi
            add ecx, 4
            add edx, 4
            jmp loop_r
        end_r:
            addss xmm0, xmm5; sommo energie ciclo q e ciclo r
            
            ;restituisco il risultato (float) sullo stack
            sub esp, 4
            movss [esp], xmm6
            fld dword [esp]
            add esp, 4
        
        pop edi
        pop esi
        mov esp, ebp
        pop ebp

        ret

            
            

