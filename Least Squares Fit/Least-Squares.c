#include"matrixandvectors.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

matrix* Least_Square_fit(vector* x, vector* y, vector* dy, int numfunc, double func(int n, double x), vector*c){
assert((*x).size==(*y).size && (*x).size==(*dy).size && numfunc>0 && (*c).size==numfunc);

int N = (*x).size;
matrix* A = matrix_alloc(N,numfunc);
matrix* R = matrix_alloc(numfunc,numfunc);
matrix* invR = matrix_alloc(numfunc,numfunc);
vector* b = vector_alloc(N);
for(int i=0; i<N;i++){
	vector_set(b,i,vector_get(y,i)/vector_get(dy,i));
}
for(int i=0; i<N;i++){
double xi = vector_get(x,i);
double dyi = vector_get(dy,i);
for(int j=0;j<numfunc;j++){
	matrix_set(A,i,j,func(j,xi)/dyi);
}}

Gram_Schmidt_QR_Decomposition(A,R);
Gram_Schmidt_QR_Solve(A,R,b,c);
Gram_Schmidt_QR_Inverse(A,R,invR);
invR = matrix_multiplication(invR,' ',invR,'T');

matrix_free(A);
matrix_free(R);
vector_free(b);
return invR;
}

void Singular_Value_Decomposition(matrix*A,matrix* S,matrix* V){
assert((*A).size1>(*A).size2 && (*S).size1==(*A).size2 && (*S).size2==(*A).size2 && (*V).size1==(*A).size2 && (*V).size2==(*A).size2);
vector* x = vector_alloc((*A).size2);
matrix* ACOPY =matrix_alloc((*A).size1,(*A).size2);
matrix* ATA =matrix_alloc((*A).size2,(*A).size2);
matrix* Q =matrix_alloc((*A).size2,(*A).size2);
matrix* D = matrix_alloc((*x).size,(*x).size);


Q=matrix_multiplication(A,'T',A,' ');
ATA=matrix_multiplication(A,'T',A,' ');

Jacobi_Diagonalisation(Q,x,V);

D = matrix_multiplication(V,'T',ATA,' ');
D = matrix_multiplication(D,' ',V,' ');

matrix_identity(S);
for(int i=0;i<(*S).size1;i++){
double d = matrix_get(D,i,i);
matrix_set(S,i,i,sqrt(d));
}

for(int i=0;i<(*D).size1;i++){
double dm = matrix_get(D,i,i);
matrix_set(D,i,i,1/sqrt(dm));
}

ACOPY = matrix_multiplication(A,' ',V,' ');
ACOPY = matrix_multiplication(A,' ', D,' ');

for(int i=0;i<(*A).size1;i++){
for(int j=0;j<(*A).size2;j++){
matrix_set(A,i,j,matrix_get(ACOPY,i,j));
}}

matrix_free(ATA);
matrix_free(Q);
matrix_free(D);
vector_free(x);
}

matrix* Least_Square_fit_SVD(vector* x, vector* y, vector* dy, int numfunc, double func(int n, double x), vector* c){

int N = (*x).size;
matrix* A = matrix_alloc(N,numfunc);
matrix* ATA = matrix_alloc(numfunc,numfunc);
matrix* U = matrix_alloc(N,numfunc);

matrix* S = matrix_alloc(numfunc,numfunc);
matrix* Sigma = matrix_alloc(numfunc,numfunc);
matrix* invS = matrix_alloc(numfunc,numfunc);
matrix* V = matrix_alloc(numfunc,numfunc);

vector* b = vector_alloc(N);
vector* cc = vector_alloc(numfunc);

for(int i=0; i<N;i++){
	vector_set(b,i,vector_get(y,i)/vector_get(dy,i));
}

for(int i=0; i<N;i++){
double xi = vector_get(x,i);
double dyi = vector_get(dy,i);
for(int j=0;j<numfunc;j++){
	matrix_set(A,i,j,func(j,xi)/dyi);
}}

ATA = matrix_multiplication(A,'T',A,' ');

Jacobi_Diagonalisation(ATA,x,V);

for(int i=0; i<numfunc;i++){
	matrix_set(S,i,i,sqrt(vector_get(x,i)));
	matrix_set(invS,i,i,1/sqrt(vector_get(x,i)));
	matrix_set(Sigma,i,i,1/(vector_get(x,i)));
}

U = matrix_multiplication(V,' ',invS,' ');
U = matrix_multiplication(A,' ',U,' ');

Gram_Schmidt_QR_Solve(U,S,b,c);
cc = matrix_vector_multiplication(V,' ',c);
for(int i=0;i<(*c).size;i++)vector_set(c,i,vector_get(cc,i));


Sigma = matrix_multiplication(V,' ',Sigma,' ');
Sigma = matrix_multiplication(Sigma,' ',V,'T');


matrix_free(A);
matrix_free(ATA);
matrix_free(invS);
matrix_free(V);
vector_free(b);
vector_free(cc);
return Sigma;
}




