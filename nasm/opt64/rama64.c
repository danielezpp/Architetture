#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <xmmintrin.h>

// Dichiarazione della funzione Assembly
extern double rama_energy64(double* phi, double* psi, int q, int r);

double min(double a, double b) {
	return (a < b) ? a : b;
}

double rama_energy(double* phi, double* psi, int n)
{
    double alpha_phi = -57.8;
    double alpha_psi = -47.0;
    double beta_phi = -119.0;
    double beta_psi = 113.0;
    double energy = 0.0;
	int i;
    
	for (i = 0; i < n; i++) {
        double alpha_dist = sqrt(pow(phi[i] - alpha_phi, 2) + pow(psi[i] - alpha_psi, 2));
        double beta_dist = sqrt(pow(phi[i] - beta_phi, 2) + pow(psi[i] - beta_psi, 2));
        energy += 0.5 * min(alpha_dist, beta_dist);
    }
	return energy;
}

int main() {
    // Dimensione dell'array (n = 256, deve essere multiplo di 4)
    const int n = 259;

    // Allocazione degli array `phi` e `psi` allineati a 16 byte
    double* phi = _mm_malloc(n*sizeof(double), 32);
    double* psi = _mm_malloc(n*sizeof(double), 32);

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
        phi[i] = -180.0f + (double)(rand() % 360); // Valori casuali tra -180 e 180
        psi[i] = -180.0f + (double)(rand() % 360); // Valori casuali tra -180 e 180
    }

    int q = (int) n/4;
    int r = n%4;
    printf("Quoziente: %d, Resto: %d\n", q,r);

    // Chiamata alla funzione Assembly
    double energy = rama_energy64(phi, psi, q, r);
    double energy_c = rama_energy(phi, psi, n);

    // Stampa del risultato
    printf("L'energia calcolata in assembly è: %f\n", energy);
    printf("L'energia calcolata in C è: %f\n", energy_c);

    // Libera la memoria allocata
    _mm_free(phi);
    _mm_free(psi);

    return EXIT_SUCCESS;
}
