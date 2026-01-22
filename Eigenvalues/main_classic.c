#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <time.h>
#include"matrixandvectors.h"

int main(int argc, char** argv){

int DIM=atoi(argv[1]);
int n=DIM, m=DIM;

time_t t;
srand((unsigned) time(&t));

		matrix* A = matrix_alloc(n,m);

		matrix* E = matrix_alloc(m,m);

		vector* x = vector_alloc(n);

for(int i=0;i<(*A).size1;i++){
for(int j=0;j<(*A).size2;j++){

		double RND = (double)rand()/RAND_MAX;
		matrix_set(A,i,j,RND);
		matrix_set(A,j,i,RND);
}}

	
int sweeps = Jacobi_Diagonalisation_Classical(A,x,E);

printf("Number of Sweeps=%i\n",sweeps);

matrix_free(A);
matrix_free(E);
vector_free(x);
return 0;
}
