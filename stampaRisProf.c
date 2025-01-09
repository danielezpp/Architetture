#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16); 
}

float* alloc_char_matrix(int rows, int cols) {
	return (float*) get_block(sizeof(float),rows*cols);
}

float* load_seq(char* filename, int *n) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	
	float* data = alloc_char_matrix(rows,cols);
	status = fread(data, sizeof(float), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	
	return data;
}


int main(int argc, char** argv) {

    int n1,n2,n3;
    float* seq1 = load_seq("phi_float.ds2", &n1);
    float* seq2 = load_seq("psi_float.ds2", &n2);
	float* en = load_seq("energy_float.ds2", &n3);

	int i;
	for(i=0; i<n3; i++)
		printf("Energia: %f,", en[i]);


    printf("phi_256: [");
    for(i=0; i<n1; i++)
        printf("%f,", seq1[i]);
    
    printf("]\n");
    printf("psi_256: [");
    for(i=0; i<n2; i++)
        printf("%f,", seq2[i]);
    
    printf("]\n");

}