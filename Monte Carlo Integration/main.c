#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"matrixandvectors.h"

double func(vector* x){
	return 1.0/(1-cos(vector_get(x,0))*cos(vector_get(x,1))*cos(vector_get(x,2)))/(M_PI*M_PI*M_PI);
}

double func_square(vector* x){
	return vector_get(x,0)*vector_get(x,0)+vector_get(x,1)*vector_get(x,1)+vector_get(x,2)*vector_get(x,2);
}
double func_sin(double x, double y){
	return sin(x)*sin(y);
}
double c(double x){
	return x*x;
}
double d(double x){
	return exp(x);
}


int main(){
printf("---------------------------------------------------------------------\n");
printf("Test of Plain Monte Carlo Multi-Dimensional Integration\n");
printf("---------------------------------------------------------------------\n");
	
	srand((unsigned) time(NULL));

	vector* a = vector_alloc(3);
	vector* b = vector_alloc(3);
	double err; 
	double res;
	double dim, N;



	//Test of x*x+y*y+z*z from 0 to 1
	vector_set(a,0,0);
	vector_set(a,1,0);
	vector_set(a,2,0);
	vector_set(b,0,1);
	vector_set(b,1,1);
	vector_set(b,2,1);

	
	dim = 3;
	N = 1e+7;
	
	Monte_Carlo_Plain(func_square, a, b, dim, N,&res,&err); 
	printf("\nCalculated Integral using Plain Monte Carlo 0f x*x+y*y+z*z from x,y,z=0 to x,y,z=1:\n\n Result=%g\nReal Result= %g\nEstimated Error= %g\nReal Error= %g\n", res, 1., err, fabs(res-1.));
	



	// Test of function in assignment
	vector_set(a,0,0);
	vector_set(a,1,0);
	vector_set(a,2,0);
	vector_set(b,0,M_PI);
	vector_set(b,1,M_PI);
	vector_set(b,2,M_PI);

	
	dim = 3;
	N = 1e+7;
	double real_res=pow(tgamma(1.0/4),4)/4.0/(M_PI*M_PI*M_PI);
	
	Monte_Carlo_Plain(func, a, b, dim, N,&res,&err); 
	printf("\nCalculated Integral using Plain Monte Carlo of function 1/(1-cosxcosycosx) from x,y,z=0 to x,y,z=PI:\n\n Result=%g\nReal Result= %g\nEstimated Error= %g\nReal Error: %g\n", res, real_res, err, fabs(res-real_res));
	
	//Now we create data for show 1/sqrt(N) dependenciy 
	for(int i=10;i<1000;i+=10){
			Monte_Carlo_Plain(func_square, a, b, dim,i,&res,&err); 
			fprintf(stderr, "%i %g\n",i,err);
	}

	// Now we test 2D adaptive integrator
	double aa = 0.;
	double bb = 1.;
	double acc = 0.00001;
	double eps = 0.00001;
	err = 0;
	double result = Adaptive_2D_Integrator(func_sin,c,d,aa,bb,acc,eps,&err);

	printf("\n\n---------------------------------------------------------------------\n");
	printf("Adaptive integration of 2 variable function on the area {(x,y):a<x<b,c(x)<y<d(d)}:\n");
	printf("---------------------------------------------------------------------\n");
	printf("We integrate sin(x)sin(y) from {(x,y):0<x<1,x*x<y<exp(x)}:\n\n");
	printf("Result=%g\n",result);
	printf("Error=%g\n",err);
	printf("Result from Wolfram Alpha=0.555661");


vector_free(a);
vector_free(b);
return 0;
}
