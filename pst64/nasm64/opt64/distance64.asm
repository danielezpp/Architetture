section .text
    global get_distance64

    get_distance64:
    ;v è in rdi e w è in rsi
        vmovapd ymm0, [rdi]
        vmovapd ymm1, [rsi]
        vsubpd ymm0, ymm1
        vmulpd ymm0, ymm0

        vhaddpd ymm0, ymm0
        vperm2f128 ymm2, ymm0, ymm0, 1
        vaddsd xmm0, xmm2
        vsqrtsd xmm0, xmm0
        
        ret

        