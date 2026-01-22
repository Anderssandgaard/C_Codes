#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include<time.h>
#include"matrixandvectors.h"

double func(vector* x){
	return 1.0/(1-cos(vector_get(x,0))*cos(vector_get(x,1))*cos(vector_get(x,2)))/(M_PI*M_PI*M_PI);
}

double func_square(vector* x){
	return vector_get(x,0)*vector_get(x,0)+vector_get(x,1)*vector_get(x,1)+vector_get(x,2)*vector_get(x,2);
}


int main(){

	vector* a = vector_alloc(3);
	vector* b = vector_alloc(3);
	double err; 
	double res;
	double dim, N;

printf("---------------------------------------------------------------------\n");
printf("Test of Plain Monte Carlo Multi-Dimensional Integration Without Multi Threading\n");
printf("---------------------------------------------------------------------\n");
	
	clock_t begin_time = clock(); 
	
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
	
	clock_t end_time = clock(); 
	double time_spend = (double)(end_time-begin_time)/CLOCKS_PER_SEC;
	
	printf("Time spend without Multi threading is %g seconds \n",time_spend);
	



printf("---------------------------------------------------------------------\n");
printf("Test of Plain Monte Carlo Multi-Dimensional Integration With Multi Threading\n");
printf("---------------------------------------------------------------------\n");
	
	begin_time = clock(); 

	// Test of function in assignment
	vector_set(a,0,0);
	vector_set(a,1,0);
	vector_set(a,2,0);
	vector_set(b,0,M_PI);
	vector_set(b,1,M_PI);
	vector_set(b,2,M_PI);

	
	dim = 3;
	N = 1e+7;
	real_res=pow(tgamma(1.0/4),4)/4.0/(M_PI*M_PI*M_PI);
	
	Monte_Carlo_Plain_OMP(func, a, b, dim, N,&res,&err); 

	printf("\nCalculated Integral using Plain Monte Carlo of function 1/(1-cosxcosycosx) from x,y,z=0 to x,y,z=PI:\n\n Result=%g\nReal Result= %g\nEstimated Error= %g\nReal Error: %g\n", res, real_res, err, fabs(res-real_res));
	
	end_time = clock(); 
	time_spend = (double)(end_time-begin_time)/CLOCKS_PER_SEC;
	
	printf("Time spend without Multi threading is %g seconds \n",time_spend);
	



vector_free(a);
vector_free(b);
return 0;
}
