#include <stdlib.h>
#include <math.h>
#include <pst32c.h>

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

/*
float dist(float coords[][3], int i, int j) 
{
     return sqrt(
        pow(coords[i][0] - coords[j][0], 2) + 
        pow(coords[i][1] - coords[j][1], 2) + 
        pow(coords[i][2] - coords[j][2], 2) ); 
}*/
/*
float hydrophobic_energy(char* s, float* coords)
{
    int n = sizeof(s) / sizeof(s[0]);
    float* ca_coords;
    int i;
    int j = 0;
    float energy = 0;

    for ( i = 6; i < n/3; i+=9)
    {
        ca_coords[j] = s[i];
        j++;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            float d = dist(coords, i, j);
            if (d < 10.0) {
                energy += (hydrophobicity[s[i]] * hydrophobicity[s[j]]) / d;
            }
        }
    }
    return energy;
}*/
type* backbone(char* s, type* phi, type* psi){return NULL;}

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


float rama_energy(float* phi, float* psi)
{
    int n = sizeof(phi) / sizeof(phi[0]);
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

float electrostatic_energy(char* s, char* coords) 
{
    float energy = 0;
    int n = sizeof(s) / sizeof(s[0]);

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

float packing_energy(char* s, float* coords) {
    float energy = 0;
    int n = sizeof(s) / sizeof(s[0]);

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
    float rama_e = rama_energy(phi, psi);
    float electro_e = electrostatic_energy(s, coords);
    float packing_e = packing_energy(s, coords );

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