#include <stdio.h>
#include <xmmintrin.h>

extern float get_distance(float* v, float* w);

int main(){
    float* v = _mm_malloc(4*sizeof(float), 16);
    float* w = _mm_malloc(4*sizeof(float), 16);
    
    for(int i=0; i<3; i++){
        v[i] = i*2.3;
        w[i] = i*1.6;
    }
    v[3] = 0.0;
    w[3] = 0.0;

    float dist = get_distance(v,w);
    printf("Distance between v and w: %f\n", dist);
}



/*type get_distance(type* v1, type* v2){
    type dx = v2[0] - v1[0];
    type dy = v2[1] - v1[1];
    type dz = v2[2] - v1[2];

    return sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
}*/