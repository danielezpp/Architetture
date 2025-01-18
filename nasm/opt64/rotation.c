double* rotation (double *axis, double theta){
    double* rotated_m = alloc_matrix(3,3);
	
    if(!rotated_m){
		printf("Errore nell'allocazione di rotated_M");
	}

	 // Copia locale per preservare `axis`
    double* normalized_axis = alloc_matrix(1,4);
    prod_axis_64_pad(axis, normalized_axis, 0);
	/*double scalar_prod = prod_scal(axis, axis, 3);
    //double scalar_prod = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]); 45k di energia in più ma più efficiente
    if (scalar_prod == 0.0) {
        printf("Errore: il vettore asse ha magnitudine zero.\n");
        free(rotated_m);
        return NULL;
    }

    normalized_axis [0] = axis [0] / scalar_prod;
	normalized_axis [1] = axis [1] / scalar_prod;
	normalized_axis [2] = axis [2] / scalar_prod;*/

    // Calcola i coefficienti quaternion
    double a = taylor_cos(theta / 2.0);
    double b = -1 * normalized_axis[0] * taylor_sin(theta / 2.0);
    double c = -1 * normalized_axis[1] * taylor_sin(theta / 2.0);
    double d = -1 * normalized_axis[2] * taylor_sin(theta / 2.0);

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