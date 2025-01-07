#include <math.h>
#include <stdlib.h>

//Calcola il prodotto scalare di un vettore con sé stesso (axis*axis)
int prod_scal(float* a){
    int lenght = sizeof(a) / sizeof(float);
    int i;
    int somma;

    for (i=0; i<lenght;i++){
        somma += a[i]*a[i];
    }

    return somma;
}

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
        double term = pow(x, 2 * i + 1) / factorial(2 * i + 1);
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

//FUNZIONE PER LA MATRICE DI ROTAZIONE DI BASE
float* rotation (float *axis, float theta){
    int n = 3;
    float* rotated_m;
    rotated_m = alloc_matrix(n,n);
    int i;

    for (i=0; i<n; i++){
        axis [i] = axis[i] / prod_scal(axis);
    }

    float a = taylor_cos(theta/2.0);
    float b = -axis[0] * taylor_sin(theta / 2.0);
    float c = -axis[1] * taylor_sin(theta / 2.0);
    float d = -axis[2] * taylor_sin(theta / 2.0);

    rotated_m[0] = a * a + b * b - c * c - d * d;
    rotated_m[1] = 2 * (b * c + a * d);
    rotated_m[2] = 2 * (b * d - a * c);

    rotated_m[3] = 2 * (b * c - a * d);
    rotated_m[4] = a * a - b * b + c * c - d * d;
    rotated_m[5] = 2 * (c * d + a * b);

    rotated_m[6] = 2 * (b * d + a * c);
    rotated_m[7] = 2 * (c * d - a * b);
    rotated_m[8] = a * a - b * b - c * c + d * d;

    return rotated_m;
}

//NORMALIZZAZIONE DELLA MATRICE
void normalize(float v[3]){
    float magn = sqrt(v[0]*v[0] + v[1]*v[1]+ v[2]*v[2]);
    if(magn = 0.0){
        printf("Errore!");
        return;
    }

    v[0] /= magn;
    v[1] /= magn;
    v[2] /= magn;
}

//PRODOTTO TRA VETTORE E MATRICE (CASO SPECIFICO V = {0,X,0})
float* prod_mat(float dist, float* rot){
    float newv[3];
    newv[0] = dist*rot[3];
    newv[1] = dist*rot[4];
    newv[2] = dist*rot[5];
    return newv;
}

//FUNZIONE BACKBONE

float* backbone(char *sequence, float *phi, float *psi){
    int n;
    n = sizeof(sequence);

    float coords[n*9];

    //Distanze standard nel backbone
    float r_ca_n= 1.46;
    float r_ca_c= 1.52;
    float r_c_n= 1.33;

    //Angoli standard nel backbone
    float th_ca_c_n= 2.028;
    float th_c_n_ca= 2.124;
    float th_n_ca_c= 1.940;

    coords [0] = 0;
    coords [1] = 0;
    coords [2] = 0;
    coords [3] = r_ca_n;
    coords [4] = 0;
    coords [5] = 0;

    int i;

    for (i=0;i<n;i++){
        int index = i*9;
        float* newv;
        float* rot;

        if (i>0)
        {
            //Posiziona N usando l'ultimo C
            float v1[3];

            v1[0]= coords[index-3] - coords[index-6];
            v1[1]= coords[index-2] - coords[index-5];
            v1[2]= coords[index-1] - coords[index-4];

            normalize(v1);

            rot = rotation(v1, th_c_n_ca);
            newv = prod_mat (r_c_n, rot);
            
            coords[index]= coords[index-3]+newv[0];
            coords[index+1]= coords[index-2]+newv[1];
            coords[index+2]= coords[index-1]+newv[2];

            //Posiziona Ca usando Phi
            float v2[3];
            
            v2[0]= coords[index] - coords[index-3];
            v2[1]= coords[index+1] - coords[index-2];
            v2[2]= coords[index+2] - coords[index-1];

            normalize(v2);

            rot = rotation(v2, phi[i]);
            newv = prod_mat (r_ca_n, rot);
            
            coords[index+3]= coords[index]+newv[0];
            coords[index+4]= coords[index+1]+newv[1];
            coords[index+5]= coords[index+2]+newv[2];
        }

        //Posiziona C usando Psi
        float v3[3];

        v3[0] = coords[index+3] - coords[index];
        v3[1] = coords[index+4] - coords[index+1];
        v3[2] = coords[index+5] - coords[index+2];
        
        normalize(v3);
        
        rot = rotation(v3, psi[i]);
        newv = prod_mat(r_ca_c, rot);
        
        coords[index+6]= coords[index+3]+newv[0];
        coords[index+7]= coords[index+4]+newv[1];
        coords[index+8]= coords[index+5]+newv[2];
    
        free(rot);
    }

    return coords;
}

