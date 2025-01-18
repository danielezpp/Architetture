section .text
    global prod_axis_64_pad
    global prod_axis_64

    prod_axis_64_pad:
       ; per calling convention in rdi ho *axis, in rsi *norm 
       ; e in rdx la flag
        vmovapd ymm0, [rdi] ; carico 4 valori di axis in ymm0
        vmovapd ymm1, ymm0 ; copio i valori in ymm1
        vmulpd ymm0, ymm0 ; li elevo al quadrato
        ;hadd
        vhaddpd ymm0, ymm0
        vperm2f128 ymm2, ymm0, ymm0, 1 ;faccio una permutazione del registro ymm0
                                    ; per poter sommare ulteriormente
        vaddsd xmm0, xmm2 ; completiamo la somma per ottenere il prodotto scalare
        vbroadcastsd ymm0, xmm0 ; propaghiamo il prodotto scalare in ymm0. Rallenta molto l'esecuzione!!
        ;vinsertf128 ymm0, ymm0, xmm0, 1

        cmp rdx, 0
        je end
        vsqrtpd ymm0, ymm0
        end:
            vdivpd ymm1, ymm0 ; dividiamo axis per prodotto_scalare

        vmovapd [rsi], ymm1 ; vado in memoria e metto il risultato in norm

        ret
    prod_axis_64:
        ;si potrebbe fare senza padding: da svolgere
        ret

