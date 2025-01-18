#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <xmmintrin.h>

// Dichiarazione della funzione Assembly
extern void rama_energy64(float* phi, float* psi, float* result);

int main() {
    // Dimensione dell'array (n = 256, deve essere multiplo di 4)
    const int n = 64;

    // Allocazione degli array `phi` e `psi` allineati a 16 byte
    float* phi = _mm_malloc(n*sizeof(float), 16);
    float* psi = _mm_malloc(n*sizeof(float), 16);

    if (!phi || !psi) {
        fprintf(stderr, "Errore nell'allocazione della memoria\n");
        return EXIT_FAILURE;
    }

      if (((uintptr_t)phi % 16 != 0) || ((uintptr_t)psi % 16 != 0)) {
        fprintf(stderr, "Errore: gli array non sono allineati a 16 byte\n");
        _mm_free(phi);
        _mm_free(psi);
        return EXIT_FAILURE;
    }

    // Inizializzazione degli array con valori di esempio
    for (int i = 0; i < n; i++) {
        phi[i] = -180.0f + (float)(rand() % 360); // Valori casuali tra -180 e 180
        printf("Phi[%d]: %f\n", i, phi[i]);
        psi[i] = -180.0f + (float)(rand() % 360); // Valori casuali tra -180 e 180
        printf("Psi[%d]: %f\n", i, psi[i]);
    }

    // Chiamata alla funzione Assembly
    float energy;
    rama_energy(phi, psi, &energy);

    // Stampa del risultato
    printf("L'energia calcolata è: %f\n", energy);

    // Libera la memoria allocata
    _mm_free(phi);
    _mm_free(psi);

    return EXIT_SUCCESS;
}
