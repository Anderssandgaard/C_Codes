#include"matrixandvectors.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

int Jacobi_Diagonalisation(matrix* A ,vector* x, matrix* V){

	int  Changed;
	int  numrot=0;
	int N=(*A).size1;

for(int i=0;i<N;i++){
	vector_set(x,i,matrix_get(A,i,i));
	}

	matrix_identity(V);

do{

	Changed=0;	

for(int i=0; i<N; i++){
	for(int j=i+1; j<N; j++){

		double Aii = vector_get(x,i);
		double Ajj = vector_get(x,j);
		double Aij = matrix_get(A,i,j);
		double theta = 0.5*atan2(2*Aij,Ajj-Aii);
		double c = cos(theta);
		double s = sin(theta);
		double Aii_prime = c*c*Aii-2*s*c*Aij+s*s*Ajj; 
		double Ajj_prime = s*s*Aii+2*s*c*Aij+c*c*Ajj;

if(Aii_prime!=Aii || Ajj_prime!=Ajj){

		Changed=1;
		numrot++;

		vector_set(x,i,Aii_prime);
		vector_set(x,j,Ajj_prime);
		matrix_set(A,i,j,0);

for(int p=0;p<i;p++){
	
		double Api = matrix_get(A,p,i);
		double Apj = matrix_get(A,p,j);
		matrix_set(A,p,i,c*Api-s*Apj);
		matrix_set(A,p,j,c*Apj+s*Api);
}

for(int p=i+1;p<j;p++){

		double Aip = matrix_get(A,i,p);
		double Apj = matrix_get(A,p,j);
		matrix_set(A,i,p,c*Aip-s*Apj);
		matrix_set(A,p,j,c*Apj+s*Aip);
}

for(int p=j+1;p<N;p++){
	
		double Aip = matrix_get(A,i,p);
		double Ajp = matrix_get(A,j,p);
		matrix_set(A,i,p,c*Aip-s*Ajp);
		matrix_set(A,j,p,c*Ajp+s*Aip);
}

for(int p=0;p<N;p++){
	
		double Vpi = matrix_get(V,p,i);
		double Vpj = matrix_get(V,p,j);
		matrix_set(V,p,i,c*Vpi-s*Vpj);
		matrix_set(V,p,j,c*Vpj+s*Vpi);
}}}}} 

while(Changed!=0);

return numrot;
}




int Jacobi_Diagonalisation_eigbyeig(matrix* A ,vector* x, matrix* V,int order,int num_of_eig){

	int  Changed;
	int  numrot=0;
	int  N=(*A).size1;

for(int i=0;i<N;i++){
	vector_set(x,i,matrix_get(A,i,i));
	}

	matrix_identity(V);

for(int i=0; i<num_of_eig; i++){

do{

	Changed=0;

for(int j=i+1; j<N; j++){

		double Aii = vector_get(x,i);
		double Ajj = vector_get(x,j);
		double Aij = matrix_get(A,i,j);
		double theta;
		if(order==0){
		theta = 0.5*atan2(2*Aij,Ajj-Aii);}
		if(order==1){
		theta = -0.5*atan2(2*Aij,-Ajj+Aii);}
		double c = cos(theta);
		double s = sin(theta);
		double Aii_prime = c*c*Aii-2*s*c*Aij+s*s*Ajj; 
		double Ajj_prime = s*s*Aii+2*s*c*Aij+c*c*Ajj;

if(Aii_prime!=Aii || Ajj_prime!=Ajj){

		Changed=1;
		numrot++;

		vector_set(x,i,Aii_prime);
		vector_set(x,j,Ajj_prime);
		matrix_set(A,i,j,0);

for(int p=0;p<i;p++){
	
		double Api = matrix_get(A,p,i);
		double Apj = matrix_get(A,p,j);
		matrix_set(A,p,i,c*Api-s*Apj);
		matrix_set(A,p,j,c*Apj+s*Api);
}

for(int p=i+1;p<j;p++){

		double Aip = matrix_get(A,i,p);
		double Apj = matrix_get(A,p,j);
		matrix_set(A,i,p,c*Aip-s*Apj);
		matrix_set(A,p,j,c*Apj+s*Aip);
}

for(int p=j+1;p<N;p++){
	
		double Aip = matrix_get(A,i,p);
		double Ajp = matrix_get(A,j,p);
		matrix_set(A,i,p,c*Aip-s*Ajp);
		matrix_set(A,j,p,c*Ajp+s*Aip);
}

for(int p=0;p<num_of_eig;p++){
	
		double Vpi = matrix_get(V,p,i);
		double Vpj = matrix_get(V,p,j);
		matrix_set(V,p,i,c*Vpi-s*Vpj);
		matrix_set(V,p,j,c*Vpj+s*Vpi);
}

}}}
while(Changed!=0);
} 
return numrot;
}



int MAXINDEX(matrix* A, int k,int N){
int m=k+1;
for(int i=k+2;i<N;i++){
	if(sqrt(matrix_get(A,k,i)*matrix_get(A,k,i))>sqrt(matrix_get(A,k,m)*matrix_get(A,k,m))){
		m=i;
	}
}
return m;
}


int Jacobi_Diagonalisation_Classical(matrix* A ,vector* x, matrix* V){

	int     N=(*A).size1;
	vector* index = vector_alloc(N);
	int     numrot=0;
	matrix_identity(V);
	int Changed;

for(int i=0;i<N;i++){
	vector_set(index,i,MAXINDEX(A,i,N));
	vector_set(x,i,matrix_get(A,i,i));
	}



do{

Changed=0;

	int m=0;
	for(int k=1; k<N-1; k++){
		double Akindex=matrix_get(A,k,vector_get(index,k));
		double Amindex=matrix_get(A,m,vector_get(index,m));
		if(sqrt(Akindex*Akindex)>sqrt(Amindex*Amindex)){
			m=k;
		}
	}

int k=m;
int l=vector_get(index,k);


		double All = vector_get(x,l);
		double Akk = vector_get(x,k);
		double Akl = matrix_get(A,k,l);
		double theta = 0.5*atan2(2*Akl,All-Akk);
		double c = cos(theta);
		double s = sin(theta);
		double Akk_prime = c*c*Akk-2*s*c*Akl+s*s*All; 
		double All_prime = s*s*Akk+2*s*c*Akl+c*c*All;


if(sqrt(Akl*Akl)>10e-5){
	
Changed=1;

numrot++;

matrix_set(A,k,l,0.0);
vector_set(x,k,Akk_prime);
vector_set(x,l,All_prime);


for(int i=0;i<k;i++){

		double Aik = matrix_get(A,i,k);
		double Ail = matrix_get(A,i,l);
		matrix_set(A,i,k,c*Aik-s*Ail);
		matrix_set(A,i,l,c*Ail+s*Aik);
}

for(int i=k+1;i<l;i++){

		double Aki = matrix_get(A,k,i);
		double Ail = matrix_get(A,i,l);
		matrix_set(A,k,i,c*Aki-s*Ail);
		matrix_set(A,i,l,c*Ail+s*Aki);
}

for(int i=l+1;i<N;i++){
	
		double Aki = matrix_get(A,k,i);
		double Ali = matrix_get(A,l,i);
		matrix_set(A,k,i,c*Aki-s*Ali);
		matrix_set(A,l,i,c*Ali+s*Aki);
}


for(int p=0;p<N;p++){
	
		double Vpk = matrix_get(V,p,k);
		double Vpl = matrix_get(V,p,l);
		matrix_set(V,p,k,c*Vpk-s*Vpl);
		matrix_set(V,p,l,c*Vpl+s*Vpk);
}

vector_set(index,k,MAXINDEX(A,k,N));
vector_set(index,l,MAXINDEX(A,l,N));

}}
while(Changed!=0);	
free(index);
return numrot;

}












