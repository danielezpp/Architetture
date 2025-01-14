#include <stdio.h>
#include <math.h>

// Funzione per calcolare il fattoriale
float factorial(int n) {
    float result = 1.0;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

// Funzione per calcolare sin(x) usando la serie di Taylor
float taylor_sin(float x) {
    int terms = 4;
    float result = 0.0;
    for (int i = 0; i < terms; i++) {
        // Calcola il termine corrente
        float term = pow(x, 2 * i + 1) / factorial(2 * i + 1);
        // Aggiungi o sottrai in base alla posizione
        if (i % 2 == 0) {
            result += term; // Termini dispari positivi
        } else {
            result -= term; // Termini dispari negativi
        }
    }
    return result;
}

// Funzione per calcolare cos(x) usando la serie di Taylor
float taylor_cos(float x) {
    int terms = 4;
    float result = 0.0;
    for (int i = 0; i < terms; i++) {
        // Calcola il termine corrente
        float term = pow(x, 2 * i) / factorial(2 * i);
        // Aggiungi o sottrai in base alla posizione
        if (i % 2 == 0) {
            result += term; // Termini pari positivi
        } else {
            result -= term; // Termini pari negativi
        }
    }
    return result;
}

int main() {
    float x = 6.27; // Valore di esempio per x
    printf("sin(%f) ≈ %f\n", x, taylor_sin(x));
    printf("Valore reale: %f\n", sin(x)); // Confronto con la funzione standard
    return 0;
}

