#include <stdlib.h>
#include <math.h>
#include "pst32c.h"

#define M_E 2.718281828459045
#define M_PI 3.14159265358979323846

void gen_rnd_mat(type* v, int N){
	int i;

	for(i=0; i<N; i++){
		// Campionamento del valore + scalatura
		v[i] = (random()*2 * M_PI) - M_PI;
	}
}

#define random() (((type) rand())/RAND_MAX)

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
type* rotation (float *axis, type theta){
    int n = 3;
    type* rotated_m;
    rotated_m = alloc_matrix(n,n);
    int i;

    for (i=0; i<n; i++){
        axis [i] = axis[i] / prod_scal(axis);
    }

    type a = taylor_cos(theta/2.0);
    type b = -axis[0] * taylor_sin(theta / 2.0);
    type c = -axis[1] * taylor_sin(theta / 2.0);
    type d = -axis[2] * taylor_sin(theta / 2.0);

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
void normalize(type v[3]){
    type magn = sqrt(v[0]*v[0] + v[1]*v[1]+ v[2]*v[2]);
    if(magn = 0.0){
        printf("Errore!");
        return;
    }

    v[0] /= magn;
    v[1] /= magn;
    v[2] /= magn;
}

//PRODOTTO TRA VETTORE E MATRICE (CASO SPECIFICO V = {0,X,0})
type* prod_mat(type dist, type* rot){
    type newv[3];
    newv[0] = dist*rot[3];
    newv[1] = dist*rot[4];
    newv[2] = dist*rot[5];
    return newv;
}

//FUNZIONE BACKBONE

type* backbone(char *sequence, type *phi, type *psi){
    int n;
    n = sizeof(sequence);

    type coords[n*9];

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
        type* newv;
        type* rot;

        if (i>0)
        {
            //Posiziona N usando l'ultimo C
            type v1[3];

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
            type v2[3];
            
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
        type v3[3];

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

type hydrophobic_energy(char* sequence, MATRIX coords, int n) {
    type energy = 0.0;

    // Itera su tutte le coppie di residui
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Calcola la distanza euclidea tra i C-alpha di i e j
            int ca_index_i = i * 3 + 1; // Indice per l'atomo C-alpha di i
            int ca_index_j = j * 3 + 1; // Indice per l'atomo C-alpha di j

            type dx = coords[ca_index_i * 3] - coords[ca_index_j * 3];
            type dy = coords[ca_index_i * 3 + 1] - coords[ca_index_j * 3 + 1];
            type dz = coords[ca_index_i * 3 + 2] - coords[ca_index_j * 3 + 2];
            type dist = sqrt(dx * dx + dy * dy + dz * dz);

            // Considera solo distanze inferiori a 10.0
            if (dist < 10.0) {
                type hydro_i = hydrophobicity[i];//sbagliato, da cambiare
                type hydro_j = hydrophobicity[j];
                energy += (hydro_i * hydro_j) / dist;
            }
        }
    }
    return energy;
}


float rama_energy(float* phi, float* psi, int n)
{
    float alpha_phi = -57.8;
    float alpha_psi = -47.0;
    float beta_phi = -119.0;
    float beta_psi = 113.0;
    float energy = 0;

    for (int i = 0; i < n; i++) {
        float alpha_dist = sqrt(pow(phi[i] - alpha_phi, 2) + pow(psi[i] - alpha_psi, 2));
        float beta_dist = sqrt(pow(phi[i] - beta_phi, 2) + pow(psi[i] - beta_psi, 2));
        energy += 0.5 * fmin(alpha_dist, beta_dist);
    }
}

float electrostatic_energy(char* s, char* coords, int n) 
{
    float energy = 0;

    for (int i = 1; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (charge[i] != 0 && charge[j] != 0) {
                float dist_sq = 0;
                for (int k = 0; k < 3; k++) {
                    float diff = coords[3 * i + k] - coords[3 * j + k];
                    dist_sq += diff * diff;
                }
                float dist = sqrt(dist_sq);
                if (dist < 10.0) {
                    energy += (charge[i] * charge[j]) / (dist * 4.0);
                }
            }
        }
    }

    return energy;
}

float packing_energy(char* s, float* coords, int n) {
    float energy = 0;

    for (int i = 1; i < n; i++) {
        float density = 0;

        for (int j = 0; j < n; j++) {
            if (i != j) {
                float dist_sq = 0;
                for (int k = 0; k < 3; k++) {
                    float diff = coords[3 * i + k] - coords[3 * j + k];
                    dist_sq += diff * diff;
                }
                float dist = sqrt(dist_sq);
                if (dist < 10.0) {
                    density += volume[j] / pow(dist, 3);
                }
            }
        }

        energy += pow(volume[i] - density, 2);
    }

    return energy;
}

float energy(char* s, float* phi, float* psi) {
    type* coords = backbone(s,phi,psi);
    int n = sizeof(s) / sizeof(s[0]);
    float rama_e = rama_energy(phi, psi, n);
    float electro_e = electrostatic_energy(s, coords, n);
    float packing_e = packing_energy(s, coords, n);

    float wrama = 1.0;
    float welec = 0.2;
    float wpack = 0.3;

    return wrama * rama_e + welec * electro_e + wpack * packing_e;
}

void simulated_annealing (char* s, type t0, type alpha, int k)
{
    int n = sizeof(s) / sizeof(s[0]);
    type T = t0;
    type* phi;
    type* psi;
    gen_rnd_mat(phi,n);
    gen_rnd_mat(psi,n);
    type e = energy(s,phi,psi);
    int t = 0;
    int i;
    type delta_phi;
    type delta_psi;
    type delta_energy;
    type temp_energy;
    type P;
    type r;
    while (T>0)
    {
        i = (int) random()%(n+1);
        delta_phi = random();
        delta_psi = random();
        phi[i] += delta_phi; 
        psi[i] += delta_psi;
        temp_energy = energy(s,phi,psi);
        delta_energy = temp_energy-e;
        if (delta_energy<=0)
            e=temp_energy;
        else
        {
            P = (float)pow(M_E, -delta_energy/(k*T) );
            r = random();
            if (r<P)
                e = temp_energy;
            else
            {
                phi[i] -= delta_phi;
                psi[i] -= delta_psi;
            }
        }
        t++;
        T = t0-sqrt(alpha*t);
    }
}