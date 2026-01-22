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

vector* matrix_get_row(matrix* A, int j){
				assert(j<(*A).size1);
				vector* Aj=vector_alloc((*A).size2);
				for(int i=0;i<(*A).size2;i++){
				(*Aj).data[i]=matrix_get(A,j,i);
				}
return Aj;
}


void matrix_set_column(matrix* A, int j,vector* v){
				assert(j<=(*A).size2 && (*A).size1==(*v).size);
				for (int i = 0; i < (*A).size1; i++){
				matrix_set(A,i,j,vector_get(v,i));
}}

void matrix_set_row(matrix* A, int j,vector* v){
				assert(j<=(*A).size1 && (*A).size2==(*v).size);
				for (int i = 0; i < (*A).size2; i++){
				matrix_set(A,j,i,vector_get(v,i));
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

vector* vector_matrix_multiplication(const vector* x,matrix* A,const char T){
				if(T=='T')A=matrix_transpose(A);
				assert((*A).size2==(*x).size);
				vector* res = vector_alloc((*A).size1);
				for(int k=0;k<(*A).size2;k++){
				double s=0;
				for(int i=0;i<(*A).size1;i++){
				s+=vector_get(x,i)*matrix_get(A,i,k);
				vector_set(res,k,s);
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
}}}

matrix* vector_outer_product(vector* x,vector* y){
assert((*x).size==(*y).size);
matrix* A = matrix_alloc((*x).size,(*x).size);
for(int i=0;i<(*x).size;i++){
for(int j=0;j<(*x).size;j++){
matrix_set(A,i,j,vector_get(x,i)*vector_get(y,j));
}}
return A;}

void matrix_add(matrix* A, matrix* B){
assert((*A).size1==(*B).size1 && (*A).size2==(*B).size2);
for(int i=0;i<(*A).size1;i++){
for(int j=0;j<(*A).size2;j++){
matrix_set(A,i,j,matrix_get(A,i,j)+matrix_get(A,i,j));
}}}


void matrix_scale(matrix* A,double a){
for(int i=0;i<(*A).size1;i++){
for(int j=0;j<(*A).size2;j++){
matrix_set(A,i,j,a*matrix_get(A,i,j));
}}}





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


// ROOT FINDING 

vector* Find_Roots_Newton(void func(vector* x, vector* fx),vector* x,double dx,double eps){

int N = (*x).size;
vector* numcalls = vector_alloc(2);
matrix* J = matrix_alloc(N,N);
matrix* R = matrix_alloc(N,N);
vector* z = vector_alloc(N);
vector* Dx = vector_alloc(N);
vector* fz = vector_alloc(N);
vector* df = vector_alloc(N);
vector* fx = vector_alloc(N);

int funccalls=0;

do{
funccalls++;
func(x,fx);

	for(int j=0;j<N;j++){
	funccalls++;
	vector_set(x,j,vector_get(x,j)+dx);
	func(x,df);
	vector_add(df,fx,-1);

		for(int i=0;i<N;i++){
		matrix_set(J,i,j,vector_get(df,i)/dx);
		}
	
	vector_set(numcalls,0,vector_get(numcalls,0)+1);
	vector_set(x,j,vector_get(x,j)-dx);
	}

		Gram_Schmidt_QR_Decomposition(J,R);
		Gram_Schmidt_QR_Solve(J,R,fx,Dx);

		vector_scale(Dx,-1.0);

		double lambda=2;


	do{
	
		z=vector_copy(x);
		vector_add(z,Dx,lambda);
		func(z,fz);
		funccalls++;
		lambda/=2;

	}while(sqrt(vector_dot_product(fz,fz))>(1-lambda/2)*sqrt(vector_dot_product(fx,fx)) && lambda>1.0/128);
	
	for (int i = 0; i < N; i++)
	{
		vector_set(x,i,vector_get(z,i));
		vector_set(fx,i,vector_get(fz,i));
	}
	

}while(sqrt(vector_dot_product(Dx,Dx))>dx && sqrt(vector_dot_product(fx,fx))>eps);
		
vector_set(numcalls,1,funccalls);

matrix_free(J);
matrix_free(R);
vector_free(z);
vector_free(Dx);
vector_free(df);
vector_free(fz);
vector_free(fx);

return numcalls;
}


vector* Find_Roots_Newton_Jacobian(void func(vector* x, vector* fx),void Jacobian(vector* x, matrix* J),vector* x,double dx,double eps){

int N = (*x).size;
vector* numcalls = vector_alloc(2);
matrix* J = matrix_alloc(N,N);
matrix* R = matrix_alloc(N,N);
vector* z = vector_alloc(N);
vector* Dx = vector_alloc(N);
vector* fz = vector_alloc(N);
vector* df = vector_alloc(N);
vector* fx = vector_alloc(N);

func(x,fx);


do{
		vector_set(numcalls,0,vector_get(numcalls,0)+1);

		Jacobian(x,J);

		Gram_Schmidt_QR_Decomposition(J,R);
		Gram_Schmidt_QR_Solve(J,R,fx,Dx);

		vector_scale(Dx,-1.0);

		double lambda=2;


	do{
		

		z=vector_copy(x);
		vector_add(z,Dx,lambda);
		func(z,fz);
		lambda/=2;
		vector_set(numcalls,1,vector_get(numcalls,1)+1);

	}while(sqrt(vector_dot_product(fz,fz))>(1-lambda/2)*sqrt(vector_dot_product(fx,fx)) && lambda>1.0/128);
	
	for (int i = 0; i < N; i++)
	{
		vector_set(x,i,vector_get(z,i));
		vector_set(fx,i,vector_get(fz,i));
	}
	

}while(sqrt(vector_dot_product(Dx,Dx))>dx && sqrt(vector_dot_product(fx,fx))>eps);

matrix_free(J);
matrix_free(R);
vector_free(z);
vector_free(Dx);
vector_free(df);
vector_free(fz);
vector_free(fx);

return numcalls;
}




vector* Find_Roots_Newton_Jacobian_Refined(void func(vector* x, vector* fx),void Jacobian(vector* x, matrix* J),vector* x,double dx,double eps){

int N = (*x).size;
vector* numcalls = vector_alloc(2);
matrix* J = matrix_alloc(N,N);
matrix* R = matrix_alloc(N,N);
vector* z = vector_alloc(N);
vector* Dx = vector_alloc(N);
vector* fz = vector_alloc(N);
vector* df = vector_alloc(N);
vector* fx = vector_alloc(N);

func(x,fx);


do{
		vector_set(numcalls,0,vector_get(numcalls,0)+1);

		Jacobian(x,J);

		Gram_Schmidt_QR_Decomposition(J,R);
		Gram_Schmidt_QR_Solve(J,R,fx,Dx);

		vector_scale(Dx,-1.0);

		double lambda=1;

		double g_0 = 0.5*vector_dot_product(fx,fx);
		double gprime_0 = -vector_dot_product(fx,fx);
		double g_lambda=0;
	do{
		
		double c =(g_lambda-g_0-gprime_0*lambda)/(lambda*lambda);
		g_lambda = g_0+gprime_0*lambda+c*lambda*lambda;
		z=vector_copy(x);
		vector_add(z,Dx,lambda);
		func(z,fz);
		lambda=gprime_0/(2*c);
		vector_set(numcalls,1,vector_get(numcalls,1)+1);

	}while(sqrt(vector_dot_product(fz,fz))>(1-lambda/2)*sqrt(vector_dot_product(fx,fx)) && lambda>1.0/128);
	
	for (int i = 0; i < N; i++)
	{
		vector_set(x,i,vector_get(z,i));
		vector_set(fx,i,vector_get(fz,i));
	}
	

}while(sqrt(vector_dot_product(Dx,Dx))>dx && sqrt(vector_dot_product(fx,fx))>eps);

matrix_free(J);
matrix_free(R);
vector_free(z);
vector_free(Dx);
vector_free(df);
vector_free(fz);
vector_free(fx);

return numcalls;
}


//MINIMISATION


vector* Newton_Minimisation(double func(vector* x, vector* Gradient, matrix* Hessian),vector* x,double eps){

int N = (*x).size;
vector* numcalls = vector_alloc(2);
matrix* Hessian = matrix_alloc(N,N);
matrix* InvHessian = matrix_alloc(N,N);
vector* Gradient = vector_alloc(N);
matrix* R = matrix_alloc(N,N);
vector* z = vector_alloc(N);
vector* Dx = vector_alloc(N);
vector* fx = vector_alloc(1);
vector* fz = vector_alloc(1);

vector_set(fx,0,func(x,Gradient,Hessian));

do{
		vector_set(numcalls,0,vector_get(numcalls,0)+1);


		Gram_Schmidt_QR_Decomposition(Hessian,R);
		Gram_Schmidt_QR_Inverse(Hessian,R,InvHessian);
		Dx = matrix_vector_multiplication(InvHessian,' ',Gradient);
		vector_scale(Dx,-1.0);


		double lambda=1;
		double alpha = 10e-4;

	do{
		

		z=vector_copy(x);
		vector_add(z,Dx,lambda);
		vector_set(fz,0,func(z,Gradient,Hessian));
		lambda/=2;
		vector_set(numcalls,1,vector_get(numcalls,1)+1);

	}while(vector_get(fz,0)>vector_get(fx,0)+alpha*lambda*vector_dot_product(Dx,Gradient) && lambda>1.0/128);
	
	for (int i = 0; i < N; i++)
	{
		vector_set(x,i,vector_get(z,i));
		vector_set(fx,i,vector_get(fz,i));
;

	}
	
}while(sqrt(vector_dot_product(Gradient,Gradient))>eps);

matrix_free(Hessian);
matrix_free(InvHessian);
vector_free(Gradient);
matrix_free(R);
vector_free(z);
vector_free(Dx);
vector_free(fx);
vector_free(fz);

return numcalls;
}


vector* Newton_Minimisation_Broyden(double f(vector* x, vector* Gradient), vector* x, double eps){
	vector* numcalls = vector_alloc(2);

	double dfxtdx, norm2_dfx, norm2_Dx, prod, fx, fxpldx, l, tol=1e-4, a=1e-4;
	int n=(*x).size;
	vector* s=vector_alloc(n);
	matrix* outer=matrix_alloc(n,n);
	vector* sTHinv=vector_alloc(n);
	vector* Ddfx=vector_alloc(n);
	vector* dfx=vector_alloc(n);
	vector* dfx2=vector_alloc(n);
	matrix* invH=matrix_alloc(n,n);
	matrix_identity(invH);
	int funccalls=0;
	do{
		
		vector_set(numcalls,0,vector_get(numcalls,0)+1);
		funccalls++;
		fx=f(x,dfx);
		s=matrix_vector_multiplication(invH,' ',dfx);
		vector_scale(s,-1.);
		l=1.0;
		vector_add(x,s,l);
		fxpldx=f(x,dfx2);
		funccalls++;

		dfxtdx = vector_dot_product(s,dfx);
		norm2_Dx=sqrt(vector_dot_product(s,s));
		while(fxpldx>fx+a*l*dfxtdx){
			l/=2;
			if(l*l*norm2_Dx<tol){
				matrix_identity(invH);
				break;
			}
			vector_add(x,s,(-1.0)*l);
			fxpldx=f(x,dfx2);
			funccalls++;

		}
		norm2_dfx=sqrt(vector_dot_product(dfx2,dfx2));

		vector_add(dfx2,dfx,-1.);
		Ddfx = matrix_vector_multiplication(invH,' ',dfx2);
		vector_scale(Ddfx,-1.);
		vector_add(Ddfx,s,1.);
		sTHinv = vector_matrix_multiplication(s,invH,' ');
		outer = vector_outer_product(Ddfx,sTHinv);
		s = matrix_vector_multiplication(invH,' ',s);
		prod = vector_dot_product(dfx2,s);
		matrix_scale(outer,1./prod);
		matrix_add(invH,outer);

	}while(norm2_dfx>eps*eps);

		vector_set(numcalls,1,funccalls);

	vector_free(Ddfx);
	vector_free(dfx);
	vector_free(dfx2);
	vector_free(s);
	matrix_free(invH);
	matrix_free(outer);
	vector_free(sTHinv);

	return numcalls;

}

double* alloc_vector(int cols)
{
return (double*)malloc(sizeof(double)*cols);
}

void free_vector(double* vector,int cols){
	free(vector);
}

double** alloc_matrix(int rows, int cols){
	
	double** matrix = (double**)malloc(sizeof(double*)*rows);
	for(int i=0;i<rows;i++){
		matrix[i]=alloc_vector(cols);
	}
	return matrix;
}

void free_matrix(double** matrix,int rows, int cols){
	for(int i=0;i<rows;i++){
		free_vector(matrix[i],cols);
	}
	free(matrix);
}



double** make_simplex(double x[], int d){

	double** simplex = alloc_matrix(d+1,d);
	for(int i=0;i<d+1;i++){
		for(int j=0;j<d;j++){
			simplex[i][j]=x[j];
		}
	}
	return simplex;
}

void print_simplex(double** simplex,int d){
				for (int i = 0; i < d+1; i++){
				for (int j = 0; j < d; j++){
				printf("%9.5f",simplex[i][j]);}
				printf("\n");
	}	
}

void reflection(double* highest, double* centroid, int dim, double* reflected){
for(int i=0;i<dim;i++){
	reflected[i]=2*centroid[i]-highest[i];
}}

void expansion(double* highest, double* centroid, int dim, double* expanded){
for(int i=0;i<dim;i++){
	expanded[i]=3*centroid[i]-2*highest[i];
}}

void contraction(double* highest, double* centroid, int dim, double* contracted){
for (int i = 0; i < dim; i++){
	contracted[i]=0.5*centroid[i]+0.5*highest[i];
}}

void reduction(double** simplex,int dim,int lo){
for(int k=0;k<dim+1;k++)
	if(k!=lo){
		for(int i=0;i<dim;i++){
			simplex[k][i]=0.5*(simplex[k][i]+simplex[lo][i]);
		}}
}

double distance(double* a, double* b, int dim){
	double s=0;
	for(int i=0;i<dim;i++){
		s+=pow(b[i]-a[i],2);
	}
return sqrt(s);
}

double size(double** simplex, int dim){
	double s=0;
	for(int k=1;k<dim+1;k++){
		double dist = distance(simplex[0],simplex[k],dim);
		if(dist>s) s=dist;
	}
return s;
}

void simplex_update(double** simplex, double* f_values, int d, int* hi, int* lo,double* centroid){
	*hi = 0;
	*lo = 0;
	double highest = f_values[0];
	double lowest = f_values[0];
	for(int k=1;k<d+1;k++){
		double next = f_values[k];
		if(next>highest){
			highest=next;
			*hi=k;
		}
		if(next<lowest){
			lowest=next;
			*lo=k;
		}
	}
	for(int i=0;i<d;i++){
		double s=0;
		for(int k=0;k<d+1;k++){
			if(k!=*hi){s+=simplex[k][i];}
		}
		centroid[i]=s/d;
	}
}


void simplex_initiate(double fun(double*),double** simplex, double* f_values, int d,int* hi,int* lo,double* centroid){
	for(int k=0;k<d+1;k++){
		f_values[k]=fun(simplex[k]);
	}
	simplex_update(simplex,f_values,d,hi,lo,centroid);
}

int downhill_simplex(double F(double*),double** simplex, int d, double simplex_size_goal){
	int hi;
	int lo;
	int k=0;
	double centroid[d];
	double F_value[d+1];
	double p1[d];
	double p2[d];

	simplex_initiate(F,simplex,F_value,d,&hi,&lo,centroid);
	while(size(simplex,d)>simplex_size_goal){
		simplex_update(simplex,F_value,d,&hi,&lo,centroid);
		reflection(simplex[hi],centroid,d,p1);
		double f_re=F(p1);
		if(f_re<F_value[lo]){
			expansion(simplex[hi],centroid,d,p2);
			double f_ex=F(p2);
			if(f_ex<f_re){
				for(int i=0;i<d;i++){simplex[hi][i]=p2[i];}
				F_value[hi]=f_ex;
						 }
			else{
				for(int i=0;i<d;i++){simplex[hi][i]=p1[i];}
				F_value[hi]=f_re;
						 }
		}
		else{
			if(f_re<F_value[hi]){
				for(int i=0;i<d;i++){simplex[hi][i]=p1[i];}
				F_value[hi]=f_re;
			}
			else{
				contraction(simplex[hi],centroid,d,p1);
				double f_co=F(p1);
				if(f_co<F_value[hi]){
					for(int i=0;i<d;i++){simplex[hi][i]=p1[i];}
					F_value[hi]=f_co;
				}
				else{
					reduction(simplex,d,lo);
					simplex_initiate(F,simplex,F_value,d,&hi,&lo,centroid);
				}
			}
		}
		k++;
	}
	return k;
}


// ODE DIFFERENTIAL EQUATIONS SOLVING AND INTEGRATION METHODS



void ODE_RK23_Step(double x, double h, vector* yx, void f(double x, vector* y, double* dydx),vector* yh, vector* error){

	int n=(*yx).size;
	int i;
	double k1[n],k2[n],k3[n],k4[n];
	vector* yt = vector_alloc(n);

	f(x,yx,k1);
	for(i=0;i<n;i++){
		vector_set(yt,i,vector_get(yx,i)+0.5*k1[1]*h);
	}
	f(x+0.5*h,yt,k2);
	for(i=0;i<n;i++){
		vector_set(yt,i,vector_get(yx,i)+0.75*k2[i]*h);
	}
	f(x+0.75*h,yt,k3);
	for(i=0;i<n;i++){
		vector_set(yh,i,vector_get(yx,i)+(2./9*k1[i]+1./3*k1[i]+4./9*k1[i])*h);
	}
	f(x+h,yh,k4);
	for(i=0;i<n;i++){
		vector_set(yt,i,vector_get(yx,i)+(7./24*k1[i]+1./4*k2[i]+1./3*k3[i]+1./8*k4[i])*h);
		vector_set(error,i,vector_get(yh,i)-vector_get(yt,i));
	}
vector_free(yt);
}



int ODE_Driver(void f(double x, vector* yx, double* dydx), matrix* ylist, vector* xlist, double b, double h, double acc, double eps, int step_max){
	
    int step=0;
    int n=(*ylist).size2;
	double x; 
    double err;
    double norm_y;
    double tol;
    double a=vector_get(xlist, 0);
	vector* yh=vector_alloc(n);
	vector* dy=vector_alloc(n);
	vector* y=vector_alloc(n);


	while(vector_get(xlist,step)<b){

		x=vector_get(xlist, step);
		y=matrix_get_row(ylist,step);

		if(x+h>b){
			h=b-x;
		}		

		ODE_RK23_Step(x,h,y,f,yh,dy);

		err=sqrt(vector_dot_product(dy,dy));
		norm_y=sqrt(vector_dot_product(yh,yh));
		tol=(norm_y*eps+acc)*sqrt(h*1.0/(b-a));

		if(err<tol){
			step++;
			if(step>step_max-1){
				fprintf(stderr, "\n\nMax num of steps reached.\n\n");
				break;
			}
			vector_set(xlist, step, x+h);
			for(int i=0;i<n;i++){
			matrix_set(ylist, step,i,vector_get(yh,i));
			}
		}
		if(err>0){
			h*=pow(tol/err,0.25)*0.95;
		}
		else{
			h*=2;
		}
	}
	vector_free(yh);
	vector_free(dy);
	return step+1;
}


double ODE_Integral(void f(double x, vector* yx, double* dydx),double a,vector* init, double b,int numfunc, double h, double acc, double eps, int step_max){
 
matrix* ylist=matrix_alloc(step_max, numfunc);
vector* xlist=vector_alloc(step_max);

vector_set(xlist, 0, a);

for(int i=0;i<numfunc;i++){
matrix_set(ylist,0,i,vector_get(init,i));
}

ODE_Driver(f,ylist,xlist,b,h,acc,eps,step_max);
double s=0;
for(int i=1;i<(*ylist).size1;i++){
    s+=matrix_get(ylist,i,0)*(vector_get(xlist,i)-vector_get(xlist,i-1));
}




matrix_free(ylist);
vector_free(xlist);
return s;
}


// ADAPTIVE INTEGRATION

double Adaptive_Integration24(double f(double x),double a,double b,double acc,double eps, double f2, double f3,int nrec,double *err){

assert(nrec<1e+6);

	double f1 = f(a+(b-a)/6);
	double f4 = f(a+5*(b-a)/6);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
	double q = (f1+f2+f3+f4)/4*(b-a);
	double tol = acc+eps*fabs(Q);
	*err = fabs(Q-q);

	if(*err<tol){
		return Q;
	}
	else{
		double Q1 = Adaptive_Integration24(f,a,(a+b)/2,acc/sqrt(2.),eps,f1,f2,nrec+1,err);
		double Q2 = Adaptive_Integration24(f,(a+b)/2,b,acc/sqrt(2.),eps,f3,f4,nrec+1,err);
		return Q1+Q2;
	}
}

double Adaptive_Integration(double f(double x),double a, double b, double acc, double eps, double *err){
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	int nrec = 0;
	return Adaptive_Integration24(f,a,b,acc,eps,f2,f3,nrec,err);
}

double Adaptive_Integration_Infinity(double f(double x),double a, double b, double acc, double eps,double *err){

int acheck=isinf(-a);
int bcheck=isinf(b);

	if(acheck && bcheck){
		double g(double t){
			return f(t/(1-t*t))*(1+t*t)/((1-t*t)*(1-t*t));
		};
		double a_new = -1, b_new = 1;
		double f2 = f(a_new+2*(b_new-a_new)/6);
		double f3 = f(a_new+4*(b_new-a_new)/6);
		int nrec = 0;
		return Adaptive_Integration24(g,a_new,b_new,acc,eps,f2,f3,nrec,err);
		
	}
	else if(acheck){
		double g(double t){
			return f(b-(1-t)/t)*1/(t*t);
		};
		double a_new = 0, b_new = 1;
		double f2 = f(a_new+2*(b_new-a_new)/6);
		double f3 = f(a_new+4*(b_new-a_new)/6);
		int nrec = 0;
		return Adaptive_Integration24(g,a_new,b_new,acc,eps,f2,f3,nrec,err);
		
	}

	else if(bcheck){
		double g(double t){
			return f(a+(1-t)/t)*1/(t*t);
		};
		double a_new = 0, b_new = 1;
		double f2 = f(a_new+2*(b_new-a_new)/6);
		double f3 = f(a_new+4*(b_new-a_new)/6);
		int nrec = 0;
		return Adaptive_Integration24(g,a_new,b_new,acc,eps,f2,f3,nrec,err);
		
	}

	else {
		
		double f2 = f(a+2*(b-a)/6);
		double f3 = f(a+4*(b-a)/6);
		int nrec = 0;
		return Adaptive_Integration24(f,a,b,acc,eps,f2,f3,nrec,err);
	}
}


double Clenshaw_Curtis(double f(double x),double a, double b, double acc, double eps, double *err){
int nrec=0;
double a_new=0.;
double b_new=M_PI;
double f2;
double f3;
double Func_Clenshaw_Curtis(double theta){
	double x = (a+b)/2.+(a-b)/2.*cos(theta);
	double res = f(x)*sin(theta)*(b-a)/2;
	return res;
}
f2 = Func_Clenshaw_Curtis(a_new+(b_new-a_new)/3);
f3 = Func_Clenshaw_Curtis(a_new+2*(b_new-a_new)/3);
return Adaptive_Integration24(Func_Clenshaw_Curtis,a_new,b_new,acc,eps,f2,f3,nrec,err);
}




















