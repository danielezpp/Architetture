section .data
    align 32
    alpha_phi dq -57.8, -57.8, -57.8, -57.8 
    alpha_psi dq -47.0, -47.0, -47.0, -47.0
    beta_phi  dq -119.0, -119.0, -119.0, -119.0
    beta_psi  dq 113.0, 113.0, 113.0, 113.0
    half      dq 0.5, 0.5, 0.5, 0.5

section .text
    global rama_energy64

    rama_energy64:
        ;rdi = *phi, rsi = *psi, rdx = q, rcx = r. Facciamo q operazioni allineate ed r non allineate
        vxorpd ymm0, ymm0
        vxorpd xmm5, xmm5
        dec rdx

        loop_q:
            cmp rdx, 0
            je end_q

            vmovapd ymm1, [rdi]; ymm1 = 4 elementi di phi
            vmovapd ymm3, ymm1; copiamo ymm1, ymm3
            vmovapd ymm2, [rsi]; xmm2 = 4 elementi di psi
            vmovapd ymm4, ymm2; copiamo ymm2 in ymm4

            vsubpd ymm1, [alpha_phi]
            vmulpd ymm1, ymm1
            vsubpd ymm2, [alpha_psi]
            vmulpd ymm2, ymm2
            vaddpd ymm1, ymm2
            vsqrtpd ymm1, ymm1

            vsubpd ymm3, [beta_phi]
            vmulpd ymm3, ymm3
            vsubpd ymm4, [beta_psi]
            vmulpd ymm4, ymm4
            vaddpd ymm3, ymm4
            vsqrtpd ymm3, ymm3

            vminpd ymm1, ymm3
            vmulpd ymm1, [half]

            vaddpd ymm0, ymm1

            dec rdx
            add rsi, 32
            add rdi, 32
            jmp loop_q
        end_q:
            vhaddpd ymm0, ymm0
            vperm2f128 ymm2, ymm0, ymm0, 1
            vaddsd xmm0, xmm2 ; in xmm0 abbiamo la somma totale dell'energia, a cui vanno 
            ;sommate le energie del ciclo resto
        loop_r:
            cmp rcx, 0
            je end_r

            vmovsd xmm6, [rdi]; xmm6= elemento phi [i] non allineato
            vmovsd xmm8, xmm6
            vsubsd xmm6, [alpha_phi]; sottraiamo a phi[i] alpha phi
            vmulsd xmm6, xmm6; elevo al quadrato
            
            vmovsd xmm7, [rsi]; xmm6= elemento psi [i] non allineato
            vmovsd xmm9, xmm7
            vsubsd xmm7, [alpha_psi]; sottraiamo a psi[i] alpha psi
            vmulsd xmm7, xmm7; elevo al quadrato
            
            vaddsd xmm6, xmm7; sommo le due differenze
            vsqrtsd xmm6, xmm6; alpha dist = sqrt(diff) va in xmm6

            vsubsd xmm8, [beta_phi]
            vmulsd xmm8, xmm8
            vsubsd xmm9, [beta_psi]
            vmulsd xmm9, xmm9

            vaddsd xmm8, xmm9
            vsqrtsd xmm8, xmm8

            vminsd xmm6, xmm8
            vaddsd xmm5, xmm6

            dec rcx
            add rdi, 8
            add rsi, 8
            jmp loop_r
        end_r:
            vaddsd xmm0, xmm5 ;sommo l'energia del ciclo q e del ciclo r in xmm0

            ret

          
            
        