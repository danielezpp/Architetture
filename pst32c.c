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

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)
#define M_PI 3.14159265358979323846

type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};		// hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};		// volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};		// charge

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
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente � che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16); 
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
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
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
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
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
* 	successivi byte: elementi del vettore 		--> k numero floating-point a precisione singola
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

//Calcola il prodotto scalare di un vettore con sé stesso (axis*axis)
int prod_scal(type* a, int lenght){
    int i;
    int somma = 0;

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

//FUNZIONE PER LA MATRICE DI ROTAZIONE DI BASE
type* rotation (float *axis, type theta){
    int n = 3;
    type* rotated_m;
    rotated_m = alloc_matrix(n,n);
    int i;

	//Potrebbe essere necessario effettuare normalizzazione 
    for (i=0; i<n; i++){
        axis [i] = axis[i] / prod_scal(axis,n);
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
    if(magn < 1e-6){
        printf("Errore!");
        return;
    }

    v[0] /= magn;
    v[1] /= magn;
    v[2] /= magn;
}

//PRODOTTO TRA VETTORE E MATRICE (CASO SPECIFICO V = {0,X,0})
type* prod_mat(type dist, type* rot){
    type* newv = (type*) malloc(3*sizeof(type));

    newv[0] = dist*rot[3];
    newv[1] = dist*rot[4];
    newv[2] = dist*rot[5];
    return newv;
}

//FUNZIONE BACKBONE

type* backbone(char *sequence, type *phi, type *psi){
    int n;
    n = strlen(sequence);

    type* coords = (type*) malloc(n*9*sizeof(type));

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
		free(newv);
    }

    return coords;
}


float rama_energy(float* phi, float* psi, int* n)
{
    float alpha_phi = -57.8;
    float alpha_psi = -47.0;
    float beta_phi = -119.0;
    float beta_psi = 113.0;
    float energy = 0;

    for (int i = 0; i < *n; i++) {
        float alpha_dist = sqrt(pow(phi[i] - alpha_phi, 2) + pow(psi[i] - alpha_psi, 2));
        float beta_dist = sqrt(pow(phi[i] - beta_phi, 2) + pow(psi[i] - beta_psi, 2));
        energy += 0.5 * fmin(alpha_dist, beta_dist);
    }
	return energy;
}

type hydrophobic_energy(char* sequence, MATRIX coords, int* n) {
    type energy = 0.0;

    // Itera su tutte le coppie di residui
    for (int i = 0; i < *n; i++) {
        for (int j = i + 1; j < *n; j++) {
            // Calcola la distanza euclidea tra i C-alpha di i e j
            int ca_index_i = i * 3 + 1; // Indice per l'atomo C-alpha di i
            int ca_index_j = j * 3 + 1; // Indice per l'atomo C-alpha di j

            type dx = coords[ca_index_i * 3] - coords[ca_index_j * 3];
            type dy = coords[ca_index_i * 3 + 1] - coords[ca_index_j * 3 + 1];
            type dz = coords[ca_index_i * 3 + 2] - coords[ca_index_j * 3 + 2];
            type dist = sqrt(dx * dx + dy * dy + dz * dz);

            // Considera solo distanze inferiori a 10.0
            if (dist < 10.0) {
                type hydro_i = hydrophobicity[sequence[i]-65];//sbagliato, da cambiare
				printf("Hydro [i]: %f", hydro_i);
                type hydro_j = hydrophobicity[sequence[j]-65];
				printf("Hydro [j]: %f", hydro_j);
                energy += (hydro_i * hydro_j) / dist;
            }
        }
    }
    return energy;
}

float electrostatic_energy(char* s, float* coords, int* n) 
{
    float energy = 0;

    for (int i = 1; i < *n; i++) {
        for (int j = i + 1; j < *n; j++) {
            if (charge[s[i]-65] != 0 && charge[s[j]-65] != 0) {
                float dist_sq = 0;
                for (int k = 0; k < 3; k++) {
                    float diff = coords[3 * i + k] - coords[3 * j + k];
                    dist_sq += diff * diff;
                }
                float dist = sqrt(dist_sq);
                if (dist < 10.0) {
                    energy += (charge[s[i]-65] * charge[s[j]-65]) / (dist * 4.0);
                }
            }
        }
    }

    return energy;
}

float packing_energy(char* s, float* coords, int* n) {
    float energy = 0;

    for (int i = 1; i < *n; i++) {
        float density = 0;

        for (int j = 0; j < *n; j++) {
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

float energy(char* s, float* phi, float* psi, int* n) {
    type* coords = backbone(s,phi,psi);
    float rama_e = rama_energy(phi, psi, n);
	float hydro_e = hydrophobic_energy(s, coords, n);
    float electro_e = electrostatic_energy(s, coords, n);
    float packing_e = packing_energy(s, coords, n);

    float wrama = 1.0;
	float whydro = 0.5;
    float welec = 0.2;
    float wpack = 0.3;

    return wrama * rama_e + whydro * hydro_e + welec * electro_e + wpack * packing_e;
}

void pst(params* input){
	// --------------------------------------------------------------
	// Codificare qui l'algoritmo di Predizione struttura terziaria
	// --------------------------------------------------------------
	printf("Rama energy: %f\n", rama_energy(input->phi, input->psi, &input->N));
	printf("Hydro energy: %f\n", hydrophobic_energy(input-> seq, backbone(input-> seq, input->phi,input->psi ), &input->N));
	printf("Electro energy: %f\n", electrostatic_energy(input-> seq, backbone(input-> seq, input->phi,input->psi ), &input->N));
	printf("Pack energy: %f\n", packing_energy(input-> seq, backbone(input-> seq, input->phi,input->psi ), &input->N));
	printf("Total energy: %f\n", energy(input-> seq, input->phi, input->psi, &input->N));
}

int main(int argc, char** argv) {
	char fname_phi[256];
	char fname_psi[256];
	char* seqfilename = NULL;
	clock_t t;
	float time;
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
	time = ((float)t)/CLOCKS_PER_SEC;

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

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->phi);
	dealloc_matrix(input->psi);
	free(input);

	return 0;
}
