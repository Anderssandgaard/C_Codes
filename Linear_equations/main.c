#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include <math.h>
#include"matrixandvectors.h"

int main(void){

//PREPARATION OF VECTORS AND MATRICES. INPUT IS RANDOMIZED AND DIMENSION IS NXM.

int n=4, m=4;

time_t t;
srand((unsigned) time(&t));

		matrix* A = matrix_alloc(n,m);
		matrix* R = matrix_alloc(m,m);
		vector* b = vector_alloc(n);
		vector* x = vector_alloc(m);

for(int i=0;i<(*A).size1;i++){
for(int j=0;j<(*A).size2;j++){
	matrix_set(A,i,j,(double)rand()/RAND_MAX);
	vector_set(b,i,(double)rand()/RAND_MAX);
}}


//Part one: Decompose A->QR with Gram Schmidt ortogonalisation
printf("------------------------------------------------------\n");
printf("Gram-Schmidt orthogonalization for QR-factorisation \n");
printf("------------------------------------------------------\n");

		matrix* Q = matrix_copy(A);

Gram_Schmidt_QR_Decomposition(Q,R);

		matrix* QR=matrix_multiplication(Q,' ',R,' ');
		matrix* QTQ=matrix_multiplication(Q,'T',Q,' ');

	matrix_print("\nRandom matrix A =\n",A);
	matrix_print("\nQ(After Gram-Schmidt QR Decomposition)=\n",Q);
	matrix_print("\nR=(Right-triangular?)\n",R);
	matrix_print("\nQ*R= (Equal to A?\n)",QR);
	matrix_print("\nQ^T*Q= (Equal to identity matrix?)\n",QTQ);


//Part Two: Solve system Ax=b via QRx=b
printf("------------------------------------------------------\n");
printf("Solve System of Linear Equation via QR Decomp using Gram-Schmidt\n");
printf("------------------------------------------------------\n");

Gram_Schmidt_QR_Solve(Q,R,b,x);

		vector* Ax = matrix_vector_multiplication(QR,' ',x);

	vector_print("\nb=(random vector)\n",b);	
	vector_print("\nx=(From Gram_Schmidt_QR_Solve)\n",x);
	vector_print("\nCheck if Ax=b?\n",Ax);


//Part Three: Matrix Inversion with QR-decomp via Gram-Schmidt
printf("------------------------------------------------------\n");
printf("Matrix Inversion with QR-decomposition via Gram-Schmidt:\n");
printf("------------------------------------------------------\n");

		matrix* A_inv = matrix_alloc(m,n);

Gram_Schmidt_QR_Inverse(Q,R,A_inv);

		matrix* A_invA = matrix_multiplication(A_inv,' ',A,' ');

	matrix_print("\nInverse of A found via QR using Gram-Schmidt\n",A_inv);
	matrix_print("\nA Inverse times A =(Equal to identity matrix?)\n",A_invA);


//Part Four: Perform Givens Rotation Method for QR-decomposition
printf("------------------------------------------------------\n");
printf("QR-decomposition via Givens:\n");
printf("------------------------------------------------------\n");

		matrix* A_Givens=matrix_copy(A);
		vector* b_Givens = vector_copy(b);

Givens_QR_Decomposition(A_Givens);

	matrix_print("\nGivens QR-Decomposition of A=\n",A_Givens);

Givens_QR_Solve(A_Givens,b_Givens);

		vector* A_Givensx = matrix_vector_multiplication(A,' ',b_Givens);
	
	vector_print("\nCheck if Ax=b? (from Givens_QR_Solve)\n",A_Givensx);

//Inversion via Givens
printf("------------------------------------------------------\n");
printf("Matrix Inversion with QR-decomposition via Givens:\n");
printf("------------------------------------------------------\n");

		matrix* A_inv_Givens = matrix_alloc(m,n);

Givens_QR_Inverse(A_Givens,A_inv_Givens);

		matrix* A_invA_Givens = matrix_multiplication(A_inv_Givens,' ',A,' ');

	matrix_print("\nInverse of A found via QR using Givens\n",A_inv_Givens);
	matrix_print("\nA Inverse times A =(Equal to identity matrix?)\n",A_invA_Givens);


//Free
matrix_free(A);
matrix_free(Q);
matrix_free(R);
matrix_free(QR);
matrix_free(QTQ);
matrix_free(A_inv);
matrix_free(A_invA);
matrix_free(A_inv_Givens);
matrix_free(A_invA_Givens);
matrix_free(A_Givens);
vector_free(x);
vector_free(b);
vector_free(Ax);
vector_free(b_Givens);
return 0;
}
