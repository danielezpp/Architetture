#include <stdio.h>
#include <xmmintrin.h>

extern void prod_mat64(double dist, double* newv);

int main(){

    double* rot = _mm_malloc(9*sizeof(double), 32);
    double* newv = _mm_malloc(4*sizeof(double), 32);
    double dist = 1.5;
    int i;

    for(i=0; i<9; i++){
        rot[i]=3.0*i;
    }
    
    newv[0]= rot[3];
    newv[1]= rot[4];
    newv[2]= rot[5];
    newv[3]= 0.0;

    prod_mat64(dist, newv);

    
    for(i=0; i<3; i++){
        printf("V[%d] : %f \n", i, newv[i]);
    }
}
