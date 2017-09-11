#include <R.h>
#include <math.h>
#include <stdio.h>

void getW(double *C, int *ptN, int *ptK, double *W) {
	int N=*ptN;
	int K=*ptK;
	int i,j,k;
	int auxij,auxji,auxii,auxjj;
	double z[N][N];
	double y[N][N];
	for (i=0 ; i<N ; i++) {
		for (j=0 ; j<N ; j++) {
			z[i][j] = 0;
			y[i][j] = 0;
		}
	}
	for (i=0 ; i<N ; i++) {
		for (j=0 ; j<N ; j++) {
			for (k=0 ; k<K ; k++) {
				auxij = N*N*k+N*i+j;
				auxji = N*N*k+N*j+i;
				auxii = N*N*k+N*i+i;
				auxjj = N*N*k+N*j+j;
				z[i][j] += C[auxii]*C[auxjj];
				y[i][j] += 0.5*C[auxjj]*(C[auxij]+C[auxji]);
			}
		}
	}
	for (i=0 ; i<N-1 ; i++) {
		for (j=i+1 ; j<N ; j++) {
			auxij = N*i+j;
			auxji = N*j+i;
			W[auxij] = (z[j][i]*y[j][i] - z[i][i]*y[i][j])/(z[j][j]*z[i][i]-z[i][j]*z[i][j]);
			W[auxji] = (z[i][j]*y[i][j] - z[j][j]*y[j][i])/(z[j][j]*z[i][i]-z[i][j]*z[i][j]);
		}
	}
}

