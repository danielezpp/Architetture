#include <stdio.h>
#include <xmmintrin.h>

extern void prod_mat(float* dist, float* newv);

int main(){

    float* rot = _mm_malloc(9*sizeof(float), 16);
    float* newv = _mm_malloc(4*sizeof(float), 16);
    float dist = 1.44;
    int i;

    for(i=0; i<9; i++){
        rot[i]=1.5*i;
    }
    
    newv[0]= rot[3];
    newv[1]= rot[4];
    newv[2]=rot[5];
    newv[3]=0.0;

    prod_mat(&dist, newv);

    
    for(i=0; i<3; i++){
        printf("V[%d] : %f \n", i, newv[i]);
    }
}
