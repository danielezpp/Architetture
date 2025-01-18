#include <stdio.h>
#include <xmmintrin.h>

extern void backbone_1(float* v1, float* coords, int index);
extern void backbone_2(float* v2, float* coords, int index);
extern void backbone_3(float* v3, float* coords, int index);

int main(){

    float* coords = _mm_malloc(18*sizeof(float), 16);
    float* newv = _mm_malloc(4*sizeof(float), 16);
    int index = 9;
    int i;

    for(i=0; i<18; i++){
        coords[i]=1.5*i;
    }
    
    newv[0]= 1.0;
    newv[1]= 2.0;
    newv[2]= 1.0;
    newv[3]= 0.0;

    int offset = 6;
    backbone_3(newv, coords, index);

    
    for(i=0; i<3; i++){
        printf("V[%d] : %f \n", i, newv[i]);
    }
}

/*
section text
    global backbone_1

    backbone_1:
        push ebp
        mov ebp, esp
        push esi

        mov ecx, [ebp+8]; puntatore a coords
        mov edx, [ebp+12]; puntatore a v1
        ;mov esi, [ebp+16]; index

        movaps xmm0, [ecx];in xmm0 abbiamo x0,y0,z0,x1 +esi-12
        movaps xmm1, xmm0; +esi-24
        subps xmm1, xmm0
        movaps [edx], xmm1

        pop esi
        mov esp, ebp
        pop ebp

        ret
*/