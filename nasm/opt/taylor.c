#include <stdio.h>

extern void taylor_sin_2(float* x);  // Dichiarazione della funzione Assembly

int main() {
    float x = 2.0;  // Esempio di input
    taylor_sin_2(&x);  // Chiamata alla funzione Assembly
    printf("Risultato: %f\n", x);
    return 0;
}