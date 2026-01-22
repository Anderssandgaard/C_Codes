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
				printf("%9.3f",matrix_get(A,i,j));}
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