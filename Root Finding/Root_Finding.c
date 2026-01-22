#include"matrixandvectors.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

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
