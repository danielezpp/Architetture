#include <stdio.h>

extern float proc(float x); // Dichiarazione della funzione Assembly

int main() {
    float x __attribute__((aligned(16))) = 5.0f;
    printf("%f\n", x);
    float a = proc(x);

    printf("Res: %f \n", a);

    return 0;
}
