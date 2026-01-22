#include"matrixandvectors.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

void Gram_Schmidt_QR_Decomposition(matrix* A, matrix* R){
for(int i=0; i<(*A).size2;i++){
		double R_ii=0;
for(int k=0; k<(*A).size1;k++){
		R_ii += matrix_get(A,k,i)*matrix_get(A,k,i);
}
		R_ii=sqrt(R_ii);
		matrix_set(R,i,i,R_ii);
for(int l=0; l<(*A).size1;l++){
		matrix_set(A,l,i,matrix_get(A,l,i)/R_ii); 
}

for(int j=i+1;j<(*A).size2;j++){
		double R_ij=0;
for(int m=0;m<(*A).size1;m++){
		R_ij += matrix_get(A,m,j)*matrix_get(A,m,i);
}
		R_ij=(R_ij);
		matrix_set(R,i,j,R_ij);

for (int n = 0; n < (*A).size1; n++){
		matrix_set(A,n,j,matrix_get(A,n,j)-matrix_get(A,n,i)*R_ij);
}}}}

void BackSubstitution(matrix* R,vector* x){
		assert((*R).size1==(*R).size2 && (*R).size2==(*x).size);
for(int i=(*R).size1-1; i>=0;i--){
		double s=0;
for(int k=i+1;k<(*R).size1;k++){
		s+=matrix_get(R,i,k)*vector_get(x,k);
}
		vector_set(x,i,(vector_get(x,i)-s)/matrix_get(R,i,i));
}}

void Gram_Schmidt_QR_Solve(matrix* Q,matrix* R,vector* b,vector* x){
		vector* Rx = matrix_vector_multiplication(Q,'T',b);
		BackSubstitution(R,Rx);
for(int i=0; i<(*Rx).size;i++){
		vector_set(x,i,vector_get(Rx,i));
		}
		free(Rx);
}

void Gram_Schmidt_QR_Inverse(matrix* Q,matrix* R,matrix* B){
for(int i=0;i<(*B).size2;i++){
		vector* unit_i = vector_alloc((*B).size2);
		vector_set_unit(unit_i,i);
		vector* RA_inv_i = matrix_vector_multiplication(Q,'T',unit_i);
		BackSubstitution(R,RA_inv_i);
		matrix_set_column(B,i,RA_inv_i);
}}

void Givens_QR_Decomposition(matrix* A_Givens){
for(int q=0; q<(*A_Givens).size2; q++){
for(int p=q+1;p<(*A_Givens).size1; p++){
		double theta = atan2(matrix_get(A_Givens,p,q),matrix_get(A_Givens,q,q));
for(int k=q; k<(*A_Givens).size2; k++){
		double xq = matrix_get(A_Givens,q,k);
		double xp = matrix_get(A_Givens,p,k);
		matrix_set(A_Givens,q,k,xq*cos(theta)+xp*sin(theta));
		matrix_set(A_Givens,p,k,-xq*sin(theta)+xp*cos(theta));
}
		matrix_set(A_Givens,p,q,theta);
}}}


void Givens_QR_QTb(matrix* QR, vector* b){
for(int q=0; q<(*QR).size2;q++){
for(int p=q+1; p<(*QR).size1;p++){
		double theta = matrix_get(QR,p,q);
		double vq = vector_get(b,q);
		double vp = vector_get(b,p);
		vector_set(b,q,vq*cos(theta)+vp*sin(theta));
		vector_set(b,p,-vq*sin(theta)+vp*cos(theta));
}}}

void Givens_QR_Solve(matrix* QR, vector* b){
		Givens_QR_QTb(QR,b);
for(int i=(*QR).size2-1;i>=0;i--){
	double s=0;
for(int k=i+1;k<(*QR).size2;k++){
	s+=matrix_get(QR,i,k)*vector_get(b,k);
}
	vector_set(b,i,(vector_get(b,i)-s)/matrix_get(QR,i,i));
}
}

void Givens_QR_Inverse(matrix* QR,matrix* B){
for(int i=0;i<(*B).size2;i++){
		vector* unit_i = vector_alloc((*B).size2);
		vector_set_unit(unit_i,i);
		Givens_QR_Solve(QR,unit_i);
		vector* x_Givens=vector_alloc((*B).size1);
for (int i = 0; i < (*B).size1; ++i){
		vector_set(x_Givens,i,vector_get(unit_i,i));}
		matrix_set_column(B,i,x_Givens);
}}
