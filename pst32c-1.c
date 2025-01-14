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
// ...

// Funzione per calcolare il prodotto scalare di due vettori
type prodotto_scalare(VECTOR vettore1, VECTOR vettore2, int n) {
    type prodotto = 0.0;
    for (int i = 0; i < n; i++) {
        prodotto += vettore1[i] * vettore2[i];
    }
    return prodotto;
}

// Funzione per calcolare il fattoriale
type fattoriale(int n) {
    type risultato = 1;
    for (int i = 1; i <= n; i++) {
        risultato *= i;
    }
    return risultato;
}

// Funzione per calcolare il coseno usando la serie di Taylor
type coseno(type x) {
    return 1 - (pow(x, 2) / fattoriale(2)) 
             + (pow(x, 4) / fattoriale(4)) 
             - (pow(x, 6) / fattoriale(6));
}

// Funzione per calcolare il seno usando la serie di Taylor
type seno(type x) {
    return x - (pow(x, 3) / fattoriale(3)) 
             + (pow(x, 5) / fattoriale(5)) 
             - (pow(x, 7) / fattoriale(7));
}

// Funzione "rotation"
MATRIX rotation(VECTOR axis, type theta) {
	type ps = prodotto_scalare(axis, axis, 3);
	
	for(int i=0; i<3; i++) {
		axis[i] /= ps;
	}

	type a = coseno(theta/2.0);
	type b = -1 * axis[0] * seno(theta/2.0);
	type c = -1 * axis[1] * seno(theta/2.0);
	type d = -1 * axis[2] * seno(theta/2.0);

	MATRIX ret = alloc_matrix(3,3);
	//Prima riga
	ret[0] = pow(a,2) + pow(b,2) - pow(c,2) - pow(d,2);
	ret[1] = 2 * ((b*c)+(a*d));
	ret[2] = 2 * ((b*d)-(a*c));

	//Seconda riga
	ret[3] = 2 * ((b*c)-(a*d));
	ret[4] = pow(a,2) + pow(c,2) - pow(b,2) - pow(d,2);
	ret[5] = 2 * ((c*d)+(a*b));

	//Terza riga
	ret[6] = 2 * ((b*d)+(a*c));
	ret[7] = 2 * ((c*d)-(a*b));
	ret[8] = pow(a,2) + pow(d,2) - pow(b,2) - pow(c,2);

	return ret;
}

// Funzione che calcola il prodotto matriciale
MATRIX prodotto_matriciale(MATRIX A, MATRIX B, int rowsA, int colsA, int rowsB, int colsB) {
	if (colsA != rowsB) {
        printf("Il numero di colonne di A è diverso dal numero di righe di B!\n");
        exit(-1);
    }

    MATRIX C = alloc_matrix(rowsA, colsB);
    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) {
			type somma = 0;
            for (int k = 0; k < colsA; k++) {
                somma += A[i * colsA + k] * B[k * colsB + j];
            }
			C[i * colsB + j] = somma;
        }
    }

    return C;
}

// Funzione che calcola e posiziona le coordinate di un atomo nella matrice coords 
void posiziona_atomo(MATRIX coords, int pos_atomo, type distanza, type angolo) {
	VECTOR v = alloc_matrix(1,3);
	for(int i=0; i<3; i++) {
		v[i] = coords[pos_atomo-3+i] - coords[pos_atomo-6+i];
	}

	type ps = prodotto_scalare(v,v,3);
	type norma = sqrt(ps);

	for(int i=0; i<3; i++) {
		v[i] /= norma;
	}
	
	MATRIX rot = rotation(v,angolo);
	VECTOR tmp = alloc_matrix(1,3);
	tmp[0] = 0; tmp[1] = distanza; tmp[2] = 0;

	MATRIX newv = prodotto_matriciale(tmp,rot,1,3,3,3);
	for(int i=0; i<3; i++) {
		coords[pos_atomo+i] = coords[pos_atomo-3+i] + newv[i];
	}
}

// Funzione "backbone"
MATRIX backbone(char *s, int n, VECTOR phi, VECTOR psi) {
	// Distanze standard nel backbone (in Angstrom)
	const type r_ca_n = 1.46;   // Distanza CA-N
	const type r_ca_c = 1.52;   // Distanza CA-C
	const type r_c_n = 1.33;    // Distanza C-N

	// Angoli standard del backbone (in radianti)
	const type theta_ca_c_n = 2.028;  // Angolo CA-C-N
	const type theta_c_n_ca = 2.124;  // Angolo C-N-CA
	const type theta_n_ca_c = 1.940;  // Angolo N-CA-C

	MATRIX coords = alloc_matrix((n*3),3);

	// Posiziona il primo aminoacido:
	coords[0] = 0; 		coords[1] = 0; coords[2] = 0;		// N
	coords[3] = r_ca_n; coords[4] = 0; coords[5] = 0;		// Ca

	// Costruisce la catena
	for(int i=0; i<n; i++) {
		int idx = i * 9;		// indice base per questo aminoacido (i * 3 (atomi) * 3 (coordinate))

		if(i>0) {
			// Posiziona N usando l'ultimo C
			posiziona_atomo(coords, idx, r_c_n, theta_c_n_ca);

			// Posiziona Ca usando phi
			posiziona_atomo(coords, idx+3, r_ca_n, phi[i]);
		}

		// Posiziona C usando psi
		posiziona_atomo(coords, idx+6, r_ca_c, psi[i]);
	}

	return coords;
}

// Funzione che restituisce il minimo tra due numeri
type min(type a, type b) {
	return (a < b) ? a : b;
}

// Funzione "rama-energy"
type rama_energy(VECTOR phi, VECTOR psi, int n) {
	//Alpha
	type alpha_phi = -57.8;
	type alpha_psi = -47.0;

	//Beta
	type beta_phi = -119.0;
	type beta_psi =  113.0;

	type energy = 0.0;
	type alpha_dist;
	type beta_dist;

	for(int i=0; i<n; i++) {
		alpha_dist = sqrt(pow(phi[i]-alpha_phi,2) + pow(psi[i]-alpha_psi,2));
		beta_dist  = sqrt(pow(phi[i]-beta_phi, 2) + pow(psi[i]-beta_psi, 2));
		energy += 0.5 * min(alpha_dist,beta_dist);
	}

	return energy;
}

// Funzione che estrae dalla matrice "coords" le coordinate degli atomi Ca in una matrice di dimensione n*3
MATRIX estrai_coordinate_atomi_ca(MATRIX coords, int n) {
	MATRIX ca_coords = alloc_matrix(n,3);
	
	for(int i=0; i<n; i++) {
		int idx = i * 9;
		for(int j=0; j<3; j++) {
			ca_coords[i*3+j] = coords[idx+3+j];
		}
	}

	return ca_coords;
}

// Funzione che calcola la distanza euclidea dati due vettori formati da coordinate (x,y,z)
type distanza_euclidea(VECTOR v1, VECTOR v2) {
    type dx = v2[0] - v1[0];
    type dy = v2[1] - v1[1];
    type dz = v2[2] - v1[2];
    return sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
}

// Funzione "hydrophobic-energy"
type hydrophobic_energy(char *s, int n, MATRIX coords) {
	// Coordinate degli atomi Ca
	MATRIX ca_coords = estrai_coordinate_atomi_ca(coords, n);
	
	type energy = 0.0;
	VECTOR v1 = alloc_matrix(1,3);
	VECTOR v2 = alloc_matrix(1,3);

	for(int i=0; i<n; i++) {
		for(int j=i+1; j<n; j++) {
			
			// Creo i due vettori contenenti le coordinate dei rispettivi atomi
			for(int k=0; k<3; k++) {
				v1[k] = ca_coords[i*3+k];
				v2[k] = ca_coords[j*3+k];
			} 

			// Calcolo la distanza euclidea tra i Ca degli aminoacidi in posizione "i" e "j"
			type dist = distanza_euclidea(v1,v2);
			if(dist < 10) {
				int pos_i = s[i] - 65;
				int pos_j = s[j] - 65;
				energy += (hydrophobicity[pos_i] * hydrophobicity[pos_j]) / dist;
			}
		}
	}

	return energy;
}

// Funzione "electrostatic-energy"
type electrostatic_energy(char *s, int n, MATRIX coords) {
	// Coordinate degli atomi Ca
	MATRIX ca_coords = estrai_coordinate_atomi_ca(coords, n);

	type energy = 0.0;
	VECTOR v1 = alloc_matrix(1,3);
	VECTOR v2 = alloc_matrix(1,3);

	for(int i=0; i<n; i++) {
		for(int j=i+1; j<n; j++) {
			
			// Creo i due vettori contenenti le coordinate dei rispettivi atomi
			for(int k=0; k<3; k++) {
				v1[k] = ca_coords[i*3+k];
				v2[k] = ca_coords[j*3+k];
			} 

			type dist = distanza_euclidea(v1,v2);
			int pos_i = s[i] - 65;
			int pos_j = s[j] - 65;
				
			if(i!=j && dist<10 && charge[pos_i]!=0 && charge[pos_j]!=0) {
				energy += (charge[pos_i] * charge[pos_j]) / (dist*4.0);
			}
		}
	}

	return energy;
}

// Funzione "packing-energy"
type packing_energy(char *s, int n, MATRIX coords) {
	// Coordinate degli atomi Ca
	MATRIX ca_coords = estrai_coordinate_atomi_ca(coords, n);

	type energy = 0;
	VECTOR v1 = alloc_matrix(1,3);
	VECTOR v2 = alloc_matrix(1,3);

	for(int i=0; i<n; i++) {
		int pos_i = s[i] - 65;

		type density = 0;
		for(int j=0; j<n; j++) {
			
			if(i!=j) {
				// Creo i due vettori contenenti le coordinate dei rispettivi atomi
				for(int k=0; k<3; k++) {
					v1[k] = ca_coords[i*3+k];
					v2[k] = ca_coords[j*3+k];
				} 

				type dist = distanza_euclidea(v1,v2);
				if(dist<10) {
					int pos_j = s[j] - 65;
					density += volume[pos_j] / pow(dist,3);
				}
			}
		}
		energy += pow(volume[pos_i]-density, 2);
	}

	return energy;
}

//Funzione "energy"
type energy(char *s, int n, VECTOR phi, VECTOR psi) {
	MATRIX coords = backbone(s,n,phi,psi);

	//Calcolo delle componenti energetiche
	type rama_e = rama_energy(phi,psi,n);
	type hydro_e = hydrophobic_energy(s,n,coords);
	type elec_e = electrostatic_energy(s,n,coords);
	type pack_e = packing_energy(s,n,coords);

	//Pesi per i diversi contributi
	type w_rama = 1.0;
	type w_hydro = 0.5;
	type w_elec = 0.2;
	type w_pack = 0.3;

	//Energia totale
	type total_e = (w_rama*rama_e) + (w_hydro*hydro_e) + (w_elec*elec_e) + (w_pack*pack_e);
	return total_e;
}

void pst(params* input){
	char *s = input->seq; 		// sequenza aminoacidica
	int n = input->N;	  		// lunghezza della sequenza
    type T = input->to;	  		// temperatura iniziale
	VECTOR phi = input->phi;	// Vettore di angoli phi
	VECTOR psi = input->psi;	// Vettore di angoli psi
    type E = energy(s, n, phi, psi);

    int t = 0;
    while(T > 0) {
		
		// genera un vicino della soluzione corrente, perturbando un angolo
        int i = random() * n;	//posizione casuale tra 0 e n;

        type v_phi =  (random()* 2 * M_PI) - M_PI; //valore reale tra -pigreco e pigreco
        phi[i] += v_phi;

        type v_psi =  (random()* 2 * M_PI) - M_PI; //valore reale tra -pigreco e pigreco
        psi[i] += v_psi;

		//calcola la variazione di energia della nuova configurazione
        type current_energy = energy(s,n,phi,psi);
        type var_energy = current_energy - E;
        
		if(var_energy <= 0){
			//accetta la nuova configurazione
            E = current_energy;
        }
		else{
            //calcolo probabilità di accettazione
            type P = exp(-var_energy / (input->k * T));
            type r = random();

            if(r <= P){
				//accetta la nuova configurazione
                E = current_energy;
            }
			else{
				//rifiuta la nuova configurazione
				//ripristina phi e psi precedenti
                phi[i] -= v_phi;
                psi[i] -= v_psi;
            }
        }
		//Aggiorna la temperatura
        t++;
        T = input->to - sqrt(input->alpha * t);
    }
	
	input->e = E;
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
	// prova(input);
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
