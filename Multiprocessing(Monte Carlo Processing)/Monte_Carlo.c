#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"matrixandvectors.h"
#include<assert.h>
#include<omp.h>
#include<gsl/gsl_rng.h>

void Randx(vector* a, vector* b, vector* x, int n, gsl_rng* r){
	for(int i=0;i<n;i++){
	vector_set(x,i,vector_get(a,i)+gsl_rng_uniform(r)*(vector_get(b,i)-vector_get(a,i)));
	}	
}

void Monte_Carlo_Plain_OMP(double f(vector* x), vector* a, vector* b, int dim, int N, double *res, double *err){
double V=1;
for(int i=0;i<dim;i++){
	V*=vector_get(b,i)-vector_get(a,i);
}
double sum=0;
double sum2=0;
double sum11=0;
double sum12=0;
double sum21=0;
double sum22=0;
double fx1;
double fx2;
vector* x1=vector_alloc(dim);
vector* x2=vector_alloc(dim);

	  const gsl_rng_type *T;
	  gsl_rng *r1,*r2;
	  gsl_rng_env_setup();
	  T = gsl_rng_default;
	  r1 = gsl_rng_alloc(T);
	  r2 = gsl_rng_alloc(T);


#pragma omp parallel sections
{
	#pragma omp section
	{
			for(int i=0;i<N/2;i++){
				Randx(a,b,x1,dim,r1);
				fx1=f(x1);
				sum11+=fx1;
				sum21+=(fx1)*(fx1);
			}
	}
	#pragma omp section
	{
			for(int i=0;i<N/2;i++){
				Randx(a,b,x2,dim,r2);
				fx2=f(x2);
				sum12+=fx2;
				sum22+=fx2*fx2;
			}
	}
}
sum=sum12+(sum11);
sum2=sum21+(sum22);

double avr=sum/N;
double var = sum2/N-avr*avr;
*res = avr*V;
*err = sqrt(var/N)*V;

gsl_rng_free(r1);
gsl_rng_free(r2);
vector_free(x1);
vector_free(x2);
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

const gsl_rng_type *T;
	  gsl_rng *r;
	  gsl_rng_env_setup();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc(T);

for(int i=0;i<N;i++){
	Randx(a,b,x,dim,r);
	fx=f(x);
	sum+=fx;
	sum2+=fx*fx;
}
double avr=sum/N;
double var = sum2/N-avr*avr;
*res = avr*V;
*err = sqrt(var/N)*V;

gsl_rng_free(r);
vector_free(x);
}
