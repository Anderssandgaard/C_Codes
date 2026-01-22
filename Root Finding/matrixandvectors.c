#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"matrixandvectors.h"
#include<assert.h>

//VECTORS//
  vector* vector_alloc(int n){
  vector* v = malloc(sizeof(vector));
  (*v).size = n;
  (*v).data = malloc(n*sizeof(double));
  if( v==NULL ) fprintf(stderr,"error in vector_alloc\n");
  return v;
}

void vector_set(vector* v, int i, double vi){
		assert(i<=(*v).size && i>=0);
		(*v).data[i]=vi;
}

double vector_get(const vector* v, int i){
		assert(i<=(*v).size && i>=0);
return (*v).data[i];
}

void vector_print(const char* s, vector* v){
		printf("%s\n",s);
		for(int i=0; i<(*v).size;i++){
		printf("%9.3f \n",vector_get(v,i));
}}

double  vector_dot_product (vector* u, vector* v){
		assert((*v).size==(*u).size);
		double kii=0;
		for (int i = 0; i < (*v).size; i++) {
		double ki = vector_get(u,i)*vector_get(v,i);
		kii+=ki;}
return kii;
}

void vector_add (vector* v, vector* u, double alpha){
		assert((*v).size==(*u).size);
		for (int i = 0; i < (*v).size; i++) {
		vector_set(v,i,alpha*vector_get(u,i)+vector_get(v,i));
}} 


void vector_set_unit(vector* v,int j){
		for(int i=0;i<(*v).size;i++){
		vector_set(v,i,0);}
		vector_set(v,j,1);
}

vector* vector_copy(const vector* v){
		vector* u = vector_alloc((*v).size);
		for(int i=0;i<(*v).size;i++){
		(*u).data[i]=(*v).data[i];};
return u;
}

void vector_scale(vector* v,double r){
		for(int i=0;i<(*v).size;i++){
		vector_set(v,i,(*v).data[i]*r);
}}


void vector_free(vector* v){
		free((*v).data); free(v); }




//MATRICES//
matrix* matrix_alloc(int n, int m){
				matrix* A=(matrix*)malloc(sizeof(matrix));
				(*A).size1=n; (*A).size2=m;
				(*A).data=(double*)malloc(n*m*sizeof(double));
return A;
}


void matrix_set(matrix* A, int i, int j, double aij){ 
				assert(i>=0 && j>=0 && i<=(*A).size1 && j<=(*A).size2);
				(*A).data[i+(*A).size1*j] = aij; }

double matrix_get(matrix* A, int i, int j){ 
				double Aij =(*A).data[i+j*(*A).size1];
return Aij;
}

void matrix_print(const char* s, matrix* A){
				printf("%s\n",s);
				for (int i = 0; i < (*A).size1; i++){
				for (int j = 0; j < (*A).size2; j++){
				printf("%9.5f",matrix_get(A,i,j));}
				printf("\n");
}}

vector* matrix_get_column(matrix* A, int j){
				assert(j<(*A).size2);
				vector* Aj=vector_alloc((*A).size1);
				for(int i=0;i<(*A).size1;i++){
				(*Aj).data[i]=matrix_get(A,i,j);
				}
return Aj;
}

void matrix_set_column(matrix* A, int j,vector* v){
				assert(j<=(*A).size2 && (*A).size1==(*v).size);
				for (int i = 0; i < (*A).size1; i++){
				matrix_set(A,i,j,vector_get(v,i));
}}

matrix* matrix_transpose(matrix* A){
				matrix* ATrans=matrix_alloc((*A).size2,(*A).size1);
				for(int i=0; i<(*A).size1;i++){
				for(int j=0; j<(*A).size2;j++){
				matrix_set(ATrans,j,i,matrix_get(A,i,j));
				}}
return ATrans;
}

matrix* matrix_multiplication(matrix* A,const char TA,matrix* B,const char TB){
				if(TA=='T')A=matrix_transpose(A);
				if(TB=='T')B=matrix_transpose(B);
				assert((*A).size2==(*B).size1);
				matrix* AxB = matrix_alloc((*A).size1,(*B).size2);
				for(int i=0;i<(*AxB).size1;i++){
				for(int j=0;j<(*AxB).size2;j++){
				double s=0;
				for(int k=0;k<(*A).size2;k++){
				s+=matrix_get(A,i,k)*matrix_get(B,k,j);
				matrix_set(AxB,i,j,s);
				}}}
return AxB;
}

vector* matrix_vector_multiplication(matrix* A,const char T,const vector* x){
				if(T=='T')A=matrix_transpose(A);
				assert((*A).size2==(*x).size);
				vector* res = vector_alloc((*A).size1);
				for(int i=0;i<(*A).size1;i++){
				double s=0;
				for(int k=0;k<(*A).size2;k++){
				s+=matrix_get(A,i,k)*vector_get(x,k);
				vector_set(res,i,s);
				}}
return res;
}


matrix* matrix_copy(matrix* A){
	matrix* B=matrix_alloc((*A).size1,(*A).size2);
	for(int i=0;i<(*A).size1*(*A).size2;i++){(*B).data[i]=(*A).data[i];}
	return B;
}

void matrix_identity(matrix* A){
for(int i=0; i<(*A).size1;i++){
for(int j=0; j<(*A).size2; j++){
		if(i==j){matrix_set(A,i,j,1);}
		else{matrix_set(A,i,j,0);}
}}



}

void matrix_free(matrix* A){
				free((*A).data); free(A); }


//QR-DECOMPOSITION

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
		vector_free(Rx);
}

void Gram_Schmidt_QR_Inverse(matrix* Q,matrix* R,matrix* B){
for(int i=0;i<(*B).size2;i++){
		vector* unit_i = vector_alloc((*Q).size1);
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



// MATRIX DIAGONALISATION

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
vector_free(index);
return numrot;

}












