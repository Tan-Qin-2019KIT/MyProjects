#include "fd.h"

float	median2d(float **mat, int ny, int nx){

	int	k;
	float		*t, med=0.0;
	
	t = vector(1,nx*ny);
	
	k = nx*ny;
	memmove(&t[1], &mat[1][1], k*sizeof(float));
	quicksort(t,1,k);

	if (k%2)	med = t[k/2+1];
	else		med = (t[k/2] + t[k/2+1]) / 2.0;
	
	free_vector(t,1,nx*ny);
	return med;
}
