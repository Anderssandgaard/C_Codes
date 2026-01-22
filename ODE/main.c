#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include"matrixandvectors.h"

void diff_sin(double t, vector* y,double* dydt){
	dydt[0]=vector_get(y,1);
	dydt[1]=(-1.0)*vector_get(y,0);
}

void diff_logistic(double t, vector* y, double* dydt){
	dydt[0]=vector_get(y,0)*(1.0-vector_get(y,0));
}

int main(){

	int step_max =1000000;
	double eps=1e-6;
	double acc=1e-6;
	double stepsize = 0.001;
	double endpoint=4*M_PI;
	double x_begin=0.;



printf("---------------------------------------------------------------------\n");
printf("ODE driver on Harmonic Differential Equation with storage of data points\n");
printf("---------------------------------------------------------------------\n");

	matrix* ylist_sin=matrix_alloc(step_max, 2);
	vector* xlist_sin=vector_alloc(step_max);
		
	vector_set(xlist_sin, 0, 0.);
	matrix_set(ylist_sin,0,0,0.);
	matrix_set(ylist_sin,0,1,1.);



	int steps_sin = ODE_Driver(diff_sin,ylist_sin,xlist_sin,endpoint,stepsize,acc,eps,step_max);

	printf("Number of Steps for Solving Harmonic Diff. Eqn. is = %i\n\n",steps_sin );

	for(int i=0; i<steps_sin; i++){
		fprintf(stderr,"%g %g %g \n", vector_get(xlist_sin, i), matrix_get(ylist_sin, i, 0), sin(vector_get(xlist_sin, i)));
	}

	vector* init = vector_alloc(2);
	vector_set(init,0,0.0);
	vector_set(init,1,1.0);

	int numfunc=2;
	double integral = ODE_Integral(diff_sin,x_begin,init,endpoint,numfunc,stepsize,acc,eps,step_max);
	printf("\nValue of integral from 0 to 4PI of Harmonic Function using ODE (Should be zero) = %2.4f\n",integral);


	fprintf(stderr, "\n\n");

printf("\n\n---------------------------------------------------------------------\n");
printf("ODE driver on Logistic Differential Equation with storage of data points\n");
printf("---------------------------------------------------------------------\n");


	matrix* ylist_logistic=matrix_alloc(step_max, 2);
	vector* xlist_logistic=vector_alloc(step_max);
		
	vector_set(xlist_logistic, 0, 0.);
	matrix_set(ylist_logistic,0,0,0.5);

	int steps_logistic = ODE_Driver(diff_logistic,ylist_logistic,xlist_logistic,endpoint,stepsize,acc,eps,step_max);

	printf("Number of Steps for Solving Logistic Diff. Eqn. is = %i\n\n",steps_logistic );
	for(int i=0; i<steps_logistic; i++){
			fprintf(stderr,"%g %g %g \n", vector_get(xlist_logistic, i), matrix_get(ylist_logistic, i, 0), 1/(1+exp(-vector_get(xlist_logistic, i))));
		}

	vector* init_log = vector_alloc(1);
	vector_set(init_log,0,0.5);

	numfunc=1;
	integral = ODE_Integral(diff_logistic,x_begin,init_log,endpoint,numfunc,stepsize,acc,eps,step_max);
	printf("\nValue of integral from 0 to 4PI of Logistic Function using ODE (approx. 11.8732)) = %2.4f\n",integral);





vector_free(xlist_sin);
matrix_free(ylist_sin);
vector_free(xlist_logistic);
matrix_free(ylist_logistic);
vector_free(init);
vector_free(init_log);

return 0;
}
