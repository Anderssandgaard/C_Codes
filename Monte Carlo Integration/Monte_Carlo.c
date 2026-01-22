#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"matrixandvectors.h"
#include<assert.h>
#define RND (double)rand()/RAND_MAX

void Randx(vector* a, vector* b, vector* x, int n){

	for(int i=0;i<n;i++){
		vector_set(x,i,vector_get(a,i)+RND*(vector_get(b,i)-vector_get(a,i)));
	}	
}

void Monte_Carlo_Plain(double f(vector* x), vector* a, vector* b, int dim, int N, double *res, double *err){
double V=1;
for(int i=0;i<dim;i++){
	V*=vector_get(b,i)-vector_get(a,i);
}
double sum=0;
double sum2=0;
double fx;
vector* x=vector_alloc(dim);

for(int i=0;i<N;i++){
	Randx(a,b,x,dim);
	fx=f(x);
	sum+=fx;
	sum2+=fx*fx;
}
double avr=sum/N;
double var = sum2/N-avr*avr;
*res = avr*V;
*err = sqrt(var/N)*V;


vector_free(x);
}

double Adaptive_2D_Integrator(double f(double x, double y),double c(double x), double d(double x),double a, double b,double acc, double eps, double *err){

double a_new, b_new;

double f_outer_integral(double x){
	double f_inner_integral(double y){
		return f(x,y);
	}
	a_new = c(x); 
	b_new = d(x);
	return Adaptive_Integration_Infinity(f_inner_integral,a_new,b_new,acc,eps,err);
}
	return Adaptive_Integration_Infinity(f_outer_integral,a,b,acc,eps,err);
}