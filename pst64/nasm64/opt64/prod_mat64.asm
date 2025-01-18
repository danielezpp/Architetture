section .text 
    global prod_mat64

    prod_mat64:
        ;per Calling conventions in xmm0 abbiamo dist, mentre in rdi abbiamo *newv
        vmovapd ymm1, [rdi]
        vbroadcastsd ymm0, xmm0 ; rallenta parecchio l'esecuzione del codice!!!
        ;vinsertf128 ymm0, ymm0, xmm0, 1
        
        vmulpd ymm1, ymm0
        vmovapd [rdi], ymm1

        ret