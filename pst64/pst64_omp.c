/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2024/25
* 
* Progetto dell'algoritmo Predizione struttura terziaria proteine 221 231 a
* in linguaggio assembly x86-32 + SSE
* 
* F. Angiulli F. Fassetti S. Nisticò, novembre 2024
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf32 pst32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o pst32.o pst32c.c -o pst32c -lm && ./pst32c $pars
* 
* oppure
* 
* ./runpst32
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)
#define M_PI 3.14159265358979323846


type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};		// hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};		// volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};		// charge
type precomputed_f [] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0};

typedef struct {
	char* seq;		// sequenza di amminoacidi
	int N;			// lunghezza sequenza
	unsigned int sd; 	// seed per la generazione casuale
	type to;		// temperatura INIZIALE
	type alpha;		// tasso di raffredamento
	type k;		// costante
	VECTOR hydrophobicity; // hydrophobicity
	VECTOR volume;		// volume
	VECTOR charge;		// charge
	VECTOR phi;		// vettore angoli phi
	VECTOR psi;		// vettore angoli psi
	type e;		// energy
	int display;
	int silent;

} params;


/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (type*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (type**).
* 
* 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente � che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,32); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri typeing-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 
* 	load_seq
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*1 byte: matrix data in row-major order --> charatteri che compongono la stringa
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	
	char* data = alloc_char_matrix(rows,cols);
	status = fread(data, sizeof(char), rows*cols, fp);
	printf("Ho caricato la sequenza\n");
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o typeing-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti l'intero 1 		--> numero intero a 32 bit
* 	successivi 4 byte: numero di elementi k     --> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> k numero typeing-point a precisione singola
*/
void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	int i;
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

/*
* 	gen_rnd_mat
* 	=========
* 
*	Genera in maniera casuale numeri reali tra -pi e pi
*	per riempire una struttura dati di dimensione Nx1
* 
*/
void gen_rnd_mat(VECTOR v, int N){
	int i;

	for(i=0; i<N; i++){
		// Campionamento del valore + scalatura
		v[i] = (random()*2 * M_PI) - M_PI;
	}
}

// PROCEDURE ASSEMBLY
extern void prova(params* input);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//								FUNZIONI FATTE DA NOI!!

//NEW METHOD
type prod_scal(VECTOR v, VECTOR w, int n) {
    type prod = 0.0;
	//#pragma parallel for reduction(+:prod)
    for (int i = 0; i < n; i++) {
        prod += v[i] * w[i];
    }
    return prod;
}

// Funzione per calcolare il fattoriale
type factorial(int n) {
    type result = 1.0;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}
//FATTORIALE OTTIMIZZATO
type opt_factorial(int n){
	return precomputed_f[n];
}

type taylor_cos(type x) {
    return 1 - (pow(x, 2) / opt_factorial(2)) 
             + (pow(x, 4) / opt_factorial(4)) 
             - (pow(x, 6) / opt_factorial(6));
}

type taylor_sin(type x){
	return x - (pow(x, 3) / opt_factorial(3)) 
             + (pow(x, 5) / opt_factorial(5)) 
             - (pow(x, 7) / opt_factorial(7));
}

//FUNZIONE PER LA MATRICE DI ROTAZIONE DI BASE
type* rotation (type *axis, type theta){
    type* rotated_m = alloc_matrix(3,3);
	
    if(!rotated_m){
		printf("Errore nell'allocazione di rotated_M");
	}

	 // Copia locale per preservare `axis`
    type normalized_axis[3];
	type scalar_prod = prod_scal(axis, axis, 3);
    //type scalar_prod = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]); 45k di energia in più ma più efficiente
    if (scalar_prod == 0.0) {
        printf("Errore: il vettore asse ha magnitudine zero.\n");
        free(rotated_m);
        return NULL;
    }

    normalized_axis [0] = axis [0] / scalar_prod;
	normalized_axis [1] = axis [1] / scalar_prod;
	normalized_axis [2] = axis [2] / scalar_prod;

    // Calcola i coefficienti quaternion
    type a = taylor_cos(theta / 2.0);
    type b = -1 * normalized_axis[0] * taylor_sin(theta / 2.0);
    type c = -1 * normalized_axis[1] * taylor_sin(theta / 2.0);
    type d = -1 * normalized_axis[2] * taylor_sin(theta / 2.0);

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

type* prod_mat(type dist, type* rot) {
    if (rot == NULL) {
        fprintf(stderr, "Errore: matrice di rotazione NULL in prod_mat\n");
        exit(EXIT_FAILURE);
    }

    type* newv = alloc_matrix(1, 3);
    if (newv == NULL) {
        fprintf(stderr, "Errore: allocazione memoria fallita in prod_mat\n");
        exit(EXIT_FAILURE);
    }


    // Calcolo dei nuovi valori
    newv[0] =  dist * rot[3]; 
    newv[1] =  dist * rot[4]; 
    newv[2] =  dist * rot[5];

    return newv;
}

//NORMALIZZAZIONE DELLA MATRICE
void normalize(type v[3]){
    type magn = sqrt(v[0]*v[0] + v[1]*v[1]+ v[2]*v[2]);
    if(magn < 1e-6){
        printf("Errore!");
        return;
    }

    v[0] /= magn;
    v[1] /= magn;
    v[2] /= magn;
}

//PRODOTTO TRA VETTORE E MATRICE (CASO SPECIFICO V = {0,X,0})
void position(MATRIX coords, int index, type dist, type theta) {
	type v[3];
	
	v[0] = coords[index-3] - coords[index-6];
	v[1] = coords[index-2] - coords[index-5];
	v[2] = coords[index-1] - coords[index-4];

    normalize(v);
	
	MATRIX rot = rotation(v,theta);

	
    MATRIX newv = prod_mat(dist, rot);

    coords [index] = coords[index-3]+ newv[0];
    coords [index+1] = coords[index-2]+ newv[1];
    coords [index+2] = coords[index-1]+ newv[2];

}



//FUNZIONE BACKBONE

type* backbone(char *sequence, type *phi, type *psi, int n){

    type* coords = alloc_matrix(n, 9);

    //Distanze standard nel backbone
    type r_ca_n= 1.46;
    type r_ca_c= 1.52;
    type r_c_n= 1.33;

    //Angoli standard nel backbone
    type th_ca_c_n= 2.028;
    type th_c_n_ca= 2.124;
    type th_n_ca_c= 1.940;

    coords [0] = 0.0;
    coords [1] = 0.0;
    coords [2] = 0.0;
    coords [3] = r_ca_n;
    coords [4] = 0.0;
    coords [5] = 0.0;

    int i;

    for (i=0;i<n;i++){
        int index = i*9;
		type* newv;
		type* rot;

        if (i>0)
        {
            //Posiziona N usando l'ultimo C
            position(coords, index, r_c_n, th_c_n_ca);
            

            //Posiziona Ca usando Phi
            position(coords, index+3, r_ca_n, phi[i]);
        }

        //Posiziona C usando Psi
        position(coords, index+6, r_ca_c, psi[i]);
    }

    return coords;
}

type min(type a, type b) {
	return (a < b) ? a : b;
}

type rama_energy(type* phi, type* psi, int n)
{
    type alpha_phi = -57.8;
    type alpha_psi = -47.0;
    type beta_phi = -119.0;
    type beta_psi = 113.0;
    type energy = 0.0;


	#pragma omp parallel for reduction(+:energy)
	for (int i = 0; i < n; i++) {
        type alpha_dist = sqrt(pow(phi[i] - alpha_phi, 2) + pow(psi[i] - alpha_psi, 2));
        type beta_dist = sqrt(pow(phi[i] - beta_phi, 2) + pow(psi[i] - beta_psi, 2));
        energy += 0.5 * min(alpha_dist, beta_dist);
	}

	return energy;
}

type rama_energy_unr(type* phi, type* psi, int n)
{
    type alpha_phi = -57.8;
    type alpha_psi = -47.0;
    type beta_phi = -119.0;
    type beta_psi = 113.0;
    type energy = 0.0;
	int r = n%4;


	#pragma omp parallel for reduction(+:energy)
	for (int i = 0; i < n-r; i+=4) {
        type alpha_dist1 = sqrt(pow(phi[i] - alpha_phi, 2) + pow(psi[i] - alpha_psi, 2));
        type beta_dist1 = sqrt(pow(phi[i] - beta_phi, 2) + pow(psi[i] - beta_psi, 2));
        energy += 0.5 * min(alpha_dist1, beta_dist1);

		type alpha_dist2 = sqrt(pow(phi[i+1] - alpha_phi, 2) + pow(psi[i+1] - alpha_psi, 2));
        type beta_dist2 = sqrt(pow(phi[i+1] - beta_phi, 2) + pow(psi[i+1] - beta_psi, 2));
        energy += 0.5 * min(alpha_dist2, beta_dist2);

		type alpha_dist3 = sqrt(pow(phi[i+2] - alpha_phi, 2) + pow(psi[i+2] - alpha_psi, 2));
        type beta_dist3 = sqrt(pow(phi[i+2] - beta_phi, 2) + pow(psi[i+2] - beta_psi, 2));
        energy += 0.5 * min(alpha_dist3, beta_dist3);

		type alpha_dist4 = sqrt(pow(phi[i+3] - alpha_phi, 2) + pow(psi[i+3] - alpha_psi, 2));
        type beta_dist4 = sqrt(pow(phi[i+3] - beta_phi, 2) + pow(psi[i+3] - beta_psi, 2));
        energy += 0.5 * min(alpha_dist4, beta_dist4);
    }
	
	if(r!=0){
		for(int j =n-r;j<n;j++){
			type alpha_distr = sqrt(pow(phi[j] - alpha_phi, 2) + pow(psi[j] - alpha_psi, 2));
        	type beta_distr = sqrt(pow(phi[j] - beta_phi, 2) + pow(psi[j] - beta_psi, 2));
        	energy += 0.5 * min(alpha_distr, beta_distr);
		}
	}
	return energy;
}

MATRIX get_ca_coords(MATRIX coords, int n){
    MATRIX ca_coords = alloc_matrix(n,3);
	//#pragma parallel for
    for(int i=0; i<n; i++){
        int ca_index = i*9;
        ca_coords[i*3] = coords[ca_index+3];
        ca_coords[i*3+1] = coords[ca_index+3+1];
        ca_coords[i*3+2] = coords[ca_index+3+2];
    }

    return ca_coords;
}

type get_distance(type* v1, type* v2){
    type dx = v2[0] - v1[0];
    type dy = v2[1] - v1[1];
    type dz = v2[2] - v1[2];

    return sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
}

type hydrophobic_energy(char* sequence, type* ca_coords, int n) {
    
    type hydro_energy = 0.0;
    //MATRIX ca_coords = get_ca_coords(coords, n);
    // Itera su tutte le coppie di residui
	#pragma omp parallel for collapse (2) reduction(+:hydro_energy)
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Creo i due vettori contenenti le coordinate dei rispettivi atomi

			type v [3] = {ca_coords[i*3], ca_coords[i*3+1], ca_coords[i*3+2]};
			type w [3] = {ca_coords[j*3], ca_coords[j*3+1], ca_coords[j*3+2]};
            
			type dist = 0.0;
			dist = get_distance(v,w);
            // Considera solo distanze inferiori a 10.0
            if (dist < 10.0) {
                type hydro_i = hydrophobicity[sequence[i]-65];
                type hydro_j = hydrophobicity[sequence[j]-65];
                hydro_energy += (hydro_i * hydro_j) / dist;
            }
        }	
    }
    return hydro_energy;
}

type hydrophobic_energy_unr(char* sequence, type* ca_coords, int n) {
    
    type hydro_energy = 0.0;
	int r = n%4;
	int q = n-r;
    // Itera su tutte le coppie di residui
	#pragma omp parallel for collapse (2) reduction(+:hydro_energy)
    for (int i = 0; i < q; i+=4) {
        for (int j = i + 1; j < q; j+=4) {
            // Creo i due vettori contenenti le coordinate dei rispettivi atomi

			for(int u=0; u<4; u++){
				type v [3] = {ca_coords[(i+u)*3], ca_coords[(i+u)*3+1], ca_coords[(i+u)*3+2]};
				type w [3] = {ca_coords[(j+u)*3], ca_coords[(j+u)*3+1], ca_coords[(j+u)*3+2]};
            
				type dist = 0.0;
				dist = get_distance(v,w);
            // Considera solo distanze inferiori a 10.0
            	if (dist < 10.0) {
                	type hydro_i = hydrophobicity[sequence[i+u]-65];
                	type hydro_j = hydrophobicity[sequence[j+u]-65];
                	hydro_energy += (hydro_i * hydro_j) / dist;
            }
			}
        }	
    }
	for(int k = q; k<n;k++){
		for (int l = k+1; l<n; l++){
			type v[3]= {ca_coords[k*3], ca_coords[k*3+1], ca_coords[k*3+2]};
			type w[3]= {ca_coords[l*3], ca_coords[l*3+1], ca_coords[l*3+2]};

			type dist=0.0;
			dist = get_distance(v,w);

			if(dist<10.0){
				type hydro_k = hydrophobicity[sequence[k]-65];
                type hydro_l = hydrophobicity[sequence[l]-65];
                hydro_energy += (hydro_k * hydro_l) / dist;
			}
		}
	}
    return hydro_energy;
}

type electrostatic_energy(char* s, type* ca_coords, int n){
	type electro_energy= 0.0;
    //MATRIX ca_coords = get_ca_coords(coords, n);

	#pragma omp parallel for collapse (2) reduction(+:electro_energy)
	for(int i=0; i< n; i++){
		for(int j= i+1; j< n; j++){
			type v[3];
    		type w[3]; 
			v[0] = ca_coords[i*3];
            v[1] = ca_coords[i*3+1];
            v[2] = ca_coords[i*3+2];

            w[0] = ca_coords[j*3];
            w[1] = ca_coords[j*3+1];
            w[2] = ca_coords[j*3+2];

            type dist = 0.0;
			dist = get_distance(v,w);

			if(dist < 10.0 && charge[s[i]-65]!=0.0 && charge[s[j]-65]!=0.0){
				electro_energy += (charge[s[i]-65]*charge[s[j]-65])/(dist*4.0);
			}
		}
	}

	return electro_energy;
}

type electrostatic_energy_unr(char* s, type* ca_coords, int n){
	type electro_energy = 0.0;
	int r = n%4;
	int q = n-r;
    // Itera su tutte le coppie di residui
	#pragma omp parallel for collapse (2) reduction(+:electro_energy)
    for (int i = 0; i < q; i+=4) {
        for (int j = i + 1; j < q; j+=4) {
            // Creo i due vettori contenenti le coordinate dei rispettivi atomi

			for(int u=0; u<4; u++){
				type v [3] = {ca_coords[(i+u)*3], ca_coords[(i+u)*3+1], ca_coords[(i+u)*3+2]};
				type w [3] = {ca_coords[(j+u)*3], ca_coords[(j+u)*3+1], ca_coords[(j+u)*3+2]};
            
				type dist = 0.0;
				dist = get_distance(v,w);
            // Considera solo distanze inferiori a 10.0
            	if(dist < 10.0 && charge[s[i+u]-65]!=0.0 && charge[s[j+u]-65]!=0.0){
				electro_energy += (charge[s[i+u]-65]*charge[s[j+u]-65])/(dist*4.0);
			}
			}
        }	
    }
	for(int k = q; k<n;k++){
		for (int l = k+1; l<n; l++){
			type v[3]= {ca_coords[k*3], ca_coords[k*3+1], ca_coords[k*3+2]};
			type w[3]= {ca_coords[l*3], ca_coords[l*3+1], ca_coords[l*3+2]};

			type dist=0.0;
			dist = get_distance(v,w);

			if(dist < 10.0 && charge[s[k]-65]!=0.0 && charge[s[l]-65]!=0.0){
				electro_energy += (charge[s[k]-65]*charge[s[l]-65])/(dist*4.0);
			}
			}
		}
	return electro_energy;
	}
    


type packing_energy(char* s, type* ca_coords, int n){
	type pack_energy = 0.0;
    //MATRIX ca_coords = get_ca_coords(coords, n);
    
	#pragma omp parallel for reduction(+:pack_energy)
	for (int i=0;i<n;i++){
        int pos_i = s[i] - 65;
		type density = 0.0;
		#pragma omp parallel for reduction(+:density) 
		for(int j=0; j<n; j++){
			type v[3];
			type w[3];
            if(i!=j){
            v[0] = ca_coords[i*3];
            v[1] = ca_coords[i*3+1];
            v[2] = ca_coords[i*3+2];

            w[0] = ca_coords[j*3];
            w[1] = ca_coords[j*3+1];
            w[2] = ca_coords[j*3+2];

            type dist = 0.0;
			dist = get_distance(v,w);
            if(dist<10.0){
                int pos_j = s[j]- 65;
                density += volume[pos_j] / pow(dist, 3);
                }
            }
			//dealloc_matrix(v);
			//dealloc_matrix(w);
		}

		pack_energy+= pow(volume[s[i]-65]-density, 2);
	}
	return pack_energy;
}

//SECONDO NOI NON E' POSSIBILE UNROLLARLO
/*type packing_energy_unr(char* s, type* ca_coords, int n){
	type pack_energy = 0.0;
    //MATRIX ca_coords = get_ca_coords(coords, n);
    
	#pragma omp parallel for reduction(+:pack_energy)
	for (int i=0;i<n;i++){
        int pos_i = s[i] - 65;
		type density = 0.0;
		#pragma omp parallel for reduction(+:density) 
		for(int j=0; j<n; j++){
            if(i!=j){
				for(int u=0; u<4; u++){
					type v [3] = {ca_coords[(i+u)*3], ca_coords[(i+u)*3+1], ca_coords[(i+u)*3+2]};
					type w [3] = {ca_coords[(j+u)*3], ca_coords[(j+u)*3+1], ca_coords[(j+u)*3+2]};
            
					type dist = 0.0;
					dist = get_distance(v,w);
            		if(dist<10.0){
                	int pos_j = s[j]- 65;
                	density += volume[pos_j] / pow(dist, 3);
                	}
            	}
			//dealloc_matrix(v);
			//dealloc_matrix(w);
			}

			pack_energy+= pow(volume[s[i]-65]-density, 2);
		}
	return pack_energy;
}*/

type energy(char* s, type* phi, type* psi, int n) {
    type* coords = backbone(s,phi,psi,n);
	type* ca_coords = get_ca_coords(coords, n);

    type rama_e = rama_energy_unr(phi, psi, n);
	type hydro_e = hydrophobic_energy_unr(s, ca_coords, n);
    type electro_e = electrostatic_energy_unr(s, ca_coords, n);

    type packing_e = packing_energy(s, ca_coords, n);

    type wrama = 1.0;
	type whydro = 0.5;
    type welec = 0.2;
    type wpack = 0.3;

    return wrama * rama_e + whydro * hydro_e + welec * electro_e + wpack * packing_e;
}


void pst(params* input){
	// --------------------------------------------------------------
	// Codificare qui l'algoritmo di Predizione struttura terziaria
	// --------------------------------------------------------------
	//simulated_annealing(input-> seq, &input->to, &input->alpha, &input->k, &input->N, &input->phi, &input-> psi, &input->e);
    char* s = input->seq;
	type T = input->to;
	int n = input->N;
	VECTOR phi = input->phi;	// Vettore di angoli phi
	VECTOR psi = input->psi;
    //type* phi = (type*) malloc(*n * sizeof(type));
    //type* psi = (type*) malloc(*n * sizeof(type));
    //gen_rnd_mat(phi,*n);
    //gen_rnd_mat(psi,*n);
	
    type e = energy(s,phi,psi, n);
    int t = 0;
    int i;
    type delta_phi;
    type delta_psi;
    type delta_energy;
    type temp_energy;
    type P;
    type r;
	printf("\nSequenza :%s\n", s);
    while (T>0)
    {
        i = random()*(n);
        delta_phi = (random()*2 * M_PI) - M_PI;
        delta_psi = (random()*2 * M_PI) - M_PI;
        phi[i] += delta_phi; 
        psi[i] += delta_psi;
        temp_energy = energy(s,phi,psi,n);
        delta_energy = temp_energy-e;
        if (delta_energy<=0)
            e=temp_energy;
        else
        {
            P = exp(-delta_energy/((input-> k)*T));
            r = random();
            if (r<=P)
                e = temp_energy;
            else
            {
                phi[i] -= delta_phi;
                psi[i] -= delta_psi;
            }
        }
        t++;
        T = input->to-sqrt(input->alpha*t);
    }
	
	
	input->e = e;

}

int main(int argc, char** argv) {
	char fname_phi[256];
	char fname_psi[256];
	char* seqfilename = NULL;
	clock_t t;
	type time;
	int d;
	
	//
	// Imposta i valori di default dei parametri
	//
	params* input = malloc(sizeof(params));
	input->seq = NULL;	
	input->N = -1;			
	input->to = -1;
	input->alpha = -1;
	input->k = -1;		
	input->sd = -1;		
	input->phi = NULL;		
	input->psi = NULL;
	input->silent = 0;
	input->display = 0;
	input->e = -1;
	input->hydrophobicity = hydrophobicity;
	input->volume = volume;
	input->charge = charge;


	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//
	if(argc <= 1){
		printf("%s -seq <SEQ> -to <to> -alpha <alpha> -k <k> -sd <sd> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tSEQ: il nome del file ds2 contenente la sequenza amminoacidica\n");
		printf("\tto: parametro di temperatura\n");
		printf("\talpha: tasso di raffredamento\n");
		printf("\tk: costante\n");
		printf("\tsd: seed per la generazione casuale\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-seq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			seqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-to") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing to value!\n");
				exit(1);
			}
			input->to = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-alpha") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing alpha value!\n");
				exit(1);
			}
			input->alpha = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sd") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing seed value!\n");
				exit(1);
			}
			input->sd = atoi(argv[par]);
			par++;
		}else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(seqfilename == NULL || strlen(seqfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	input->seq = load_seq(seqfilename, &input->N, &d);
	

	
	if(d != 1){
		printf("Invalid size of sequence file, should be %ix1!\n", input->N);
		exit(1);
	} 

	if(input->to <= 0){
		printf("Invalid value of to parameter!\n");
		exit(1);
	}

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	if(input->alpha <= 0){
		printf("Invalid value of alpha parameter!\n");
		exit(1);
	}

	input->phi = alloc_matrix(input->N, 1);
	input->psi = alloc_matrix(input->N, 1);
	// Impostazione seed 
	srand(input->sd);
	// Inizializzazione dei valori
	gen_rnd_mat(input->phi, input->N);
	gen_rnd_mat(input->psi, input->N);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", seqfilename);
		printf("Sequence lenght: %d\n", input->N);
	}

	// COMMENTARE QUESTA RIGA!
	//prova(input);
	//

	//
	// Predizione struttura terziaria
	//
	t = clock();
	pst(input);
	t = clock() - t;
	time = ((type)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("PST time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname_phi, "out32_%d_%d_phi.ds2", input->N, input->sd);
	save_out(fname_phi, input->phi, input->N);
	sprintf(fname_psi, "out32_%d_%d_psi.ds2", input->N, input->sd);
	save_out(fname_psi, input->psi, input->N);
	if(input->display){
		if(input->phi == NULL || input->psi == NULL)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("energy: %f, phi: [", input->e);
			for(i=0; i<input->N; i++){
				printf("%f,", input->phi[i]);
			}
			printf("]\n");
			printf("psi: [");
			for(i=0; i<input->N; i++){
				printf("%f,", input->psi[i]);
			}
			printf("]\n");
		}
	}

	type* seq1 = load_data("phi_256_to20_k1_alpha1_sd3.ds2", &input->N, &input->N);
    type* seq2 = load_data("psi_256_to20_k1_alpha1_sd3.ds2", &input->N, &input->N);
	type* en = load_data("energy_double.ds2", &input->N, &input->N);

	int i;

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->phi);
	dealloc_matrix(input->psi);
	free(input);

	return 0;
}