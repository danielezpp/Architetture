#include <stdio.h>
#include <stdlib.h>
#include <xmmintrin.h> // Per _mm_malloc e _mm_free

extern void process_array(float* array, int n);
extern float prod_scal(float* v1, float* v2, int n);

int main() {
    int n = 8; // Dimensione dell'array
    float* v1 = (float*)_mm_malloc(n * sizeof(float), 16); // Allocazione allineata a 16 byte
    float* v2 = (float*)_mm_malloc(n * sizeof(float), 16);

    int i;
    for(i=0; i<n; i++){
        v1[i] = 1.0*i;
        v2[i] = 3.5*i;
    }

    printf("Risultato: %f\n", prod_scal(v1, v2, n));

    _mm_free(v1);
    _mm_free(v2); // Libera la memoria allineata
    return 0;
}
