#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <time.h>
#include"matrixandvectors.h"

int main(int argc, char** argv){

//PREPARATION OF VECTORS AND MATRICES. INPUT IS RANDOMIZED AND DIMENSION IS nXm.

		time_t t;
		srand((unsigned) time(&t));

		int DIM=atoi(argv[1]);
		int n=DIM, m=DIM;

		matrix* A = matrix_alloc(n,m);
		matrix* Q1 = matrix_alloc(n,m);
		matrix* Q2 = matrix_alloc(n,m);
		matrix* Q3 = matrix_alloc(n,m);
		matrix* Q4 = matrix_alloc(n,m);

		matrix* E = matrix_alloc(m,m);
		matrix* V1 = matrix_alloc(m,m);
		matrix* V2 = matrix_alloc(m,m);
		matrix* V3 = matrix_alloc(m,m);
		matrix* V4 = matrix_alloc(m,m);

		matrix* VTAV = matrix_alloc(n,m);
		matrix* VTAV2 = matrix_alloc(n,m);
		matrix* VTAV3 = matrix_alloc(n,m);
		matrix* VTAV4 = matrix_alloc(n,m);

		vector* x = vector_alloc(n);
		vector* v1 = vector_alloc(n);
		vector* v2 = vector_alloc(n);
		vector* v3 = vector_alloc(n);
		vector* v4 = vector_alloc(n);

for(int i=0;i<(*A).size1;i++){
for(int j=0;j<(*A).size2;j++){

		double RND = (double)rand()/RAND_MAX;
		matrix_set(A,i,j,RND);
		matrix_set(A,j,i,RND);
	}}

//COPY MATRICES TO PRESERVE A,E,X
		Q1=matrix_copy(A);
		Q2=matrix_copy(A);
		Q3=matrix_copy(A);
		Q4=matrix_copy(A);
		V1=matrix_copy(E);
		V2=matrix_copy(E);
		V3=matrix_copy(E);
		V4=matrix_copy(E);
		v1=vector_copy(x);
		v2=vector_copy(x);
		v3=vector_copy(x);
		v4=vector_copy(x);





printf("------------------------------------------------------\n");
printf("Jacobi Diagonalisation Using a Row by Row Algorithm \n");
printf("------------------------------------------------------\n");


int sweeps1 = Jacobi_Diagonalisation(Q1,v1,V1);

			VTAV = matrix_multiplication(A,' ',V1,' ');
			VTAV = matrix_multiplication(V1,'T',VTAV,' ');

			matrix_print("\nA(To be Diagonalised)=\n",A);
			matrix_print("\nA(After Diagonalisation)=\n",Q1);
			matrix_print("\nEigenvector Matrix for A=\n",V1);
			vector_print("\nEigenvalue Vector for A=\n",v1);
			matrix_print("\nVTAV=(Diagonal Equal to Eigenvalue Vector?)\n",VTAV);
			printf("\nNumber of Rotations= %i\n",sweeps1);


printf("\n-----------------------------------------------------\n-");
printf("Jacobi Diagonalisation Using a Eigenvalye-by-Eigenvalue Algorithm. Code is exemplified for calculation lowest eigenvalue.\n");
printf("------------------------------------------------------\n");
	
		int Number_of_eigenvalues = 1;

int sweeps2 = Jacobi_Diagonalisation_eigbyeig(Q2,v2,V2,0,Number_of_eigenvalues);

		VTAV2 = matrix_multiplication(A,' ',V2,' ');
		VTAV2 = matrix_multiplication(V2,'T',VTAV2,' ');
	
		printf("\nEigenvalues=\n");
		for(int i=0;i<Number_of_eigenvalues;i++){
		printf("%4.3f\n",vector_get(v2,i));
		}
		printf("\nNumber of Rotations= %i\n",sweeps2);



printf("\n-----------------------------------------------------\n-");
printf("Jacobi Diagonalisation Using a Eigenvalye-by-Eigenvalue Algorithm. Code is exemplified for calculation of all eigenvalues in descending order.\n");
printf("------------------------------------------------------\n");
	

int sweeps3 = Jacobi_Diagonalisation_eigbyeig(Q3,v3,V3,1,(*A).size1);

		VTAV3 = matrix_multiplication(A,' ',V3,' ');
		VTAV3 = matrix_multiplication(V3,'T',VTAV3,' ');

		matrix_print("\nEigenvector Matrix for A=\n",V3);
		vector_print("\nEigenvalue Vector for A=\n",v3);
		matrix_print("\nVTAV=(Diagonal Equal to Eigenvalue Vector?)\n",VTAV3);
		printf("\nNumber of Rotations= %i\n",sweeps3);



printf("\n-----------------------------------------------------\n-");
printf("Jacobi Diagonalisation Using a Classic Algorithm \n");
printf("------------------------------------------------------\n");
	

int sweeps4 = Jacobi_Diagonalisation_Classical(Q4,v4,V4);

		VTAV4 = matrix_multiplication(A,' ',V4,' ');
		VTAV4 = matrix_multiplication(V4,'T',VTAV4,' ');

		matrix_print("\nEigenvector Matrix for A=\n",V4);
		vector_print("\nEigenvalue Vector for A=\n",v4);
		matrix_print("\nVTAV=(Diagonal Equal to Eigenvalue Vector?)\n",VTAV4);
		printf("\nNumber of Rotations= %i\n",sweeps4);



//FREE ALLOCATIONS
matrix_free(A);
matrix_free(E);
matrix_free(V1);
matrix_free(V2);
matrix_free(V3);
matrix_free(V4);
vector_free(x);
vector_free(v1);
vector_free(v2);
vector_free(v3);
vector_free(v4);
matrix_free(VTAV);
matrix_free(VTAV2);
matrix_free(VTAV3);
matrix_free(VTAV4);
return 0;
}
