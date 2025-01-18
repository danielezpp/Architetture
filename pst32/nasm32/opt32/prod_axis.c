#include <stdio.h>
#include <xmmintrin.h>

extern void prod_axis(float* axis, float* norm);


int main(){
    float* norm = _mm_malloc(4*sizeof(float), 16);
    float* axis = _mm_malloc(4*sizeof(float), 16);
    axis [0] = 1.0;
    axis [1] = 0.8;
    axis [2] = 2.0;
    axis [3] = 0.0;
    
    prod_axis(axis, norm);
    
    int i;
    for(i=0; i<3; i++){
        printf("V[%d] : %f \n", i, norm[i]);
    }
}
/*
float normalized_axis[3];
	float scalar_prod = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    if (scalar_prod == 0.0) {
        printf("Errore: il vettore asse ha magnitudine zero.\n");
        free(rotated_m);
        return NULL;
    }

    normalized_axis [0] = axis [0] / scalar_prod;
	normalized_axis [1] = axis [1] / scalar_prod;
	normalized_axis [2] = axis [2] / scalar_prod;*/