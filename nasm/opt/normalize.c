#include <math.h>
#include <stdio.h>
#include <xmmintrin.h>

extern void normalize_v(float* v);

void normalize(float v[3]){
    float* vec = _mm_malloc(4*sizeof(float), 16);
    vec[0] = v[0];
    vec[1] = v[1];
    vec[2] = v[2];
    vec[3] = 0.0;

    normalize_v(vec);

    v[0] = vec[0];
    v[1] = vec[1];
    v[2] = vec[2];

    /*if(magn < 1e-6){
        printf("Errore!");
        return;
    }

    v[0] /= magn;
    v[1] /= magn;
    v[2] /= magn;*/
}

int main(){
    float v [] = {1.0,2.0,3.0};
    normalize(v);
    
    int i;
    for(i=0; i<3; i++){
        printf("V[%d] : %f \n", i, v[i]);
    }
}