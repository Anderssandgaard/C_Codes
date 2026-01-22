#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include"neurons.h"

//Training Functions
double train_func_log(double x){return x/(1+exp(-x));}
double derivative_train_func_log(double x){return 1/(1+exp(-x))+x/(1+exp(-x))*(1-1/(1+exp(-x)));}
double train_func_gauss(double x){return x*exp(-x*x);}
double derivative_train_func_gauss(double x){return  exp(-x*x)-2*x*x*exp(-x*x);}

//ODE to be solved
double Logistic_ODE(double x,double y){return y*(1-y);}
double Gaussian_ODE(double x,double y){return -x*y;}

int main(){

	//The span of x we are interested to solve for
	double a=-5,b=5;


////////////////////////////////////////////////////////////////////////////////
//Logistic Differential equation:
////////////////////////////////////////////////////////////////////////////////

	//Number of hidden neurons to train
	int num_neurons=20;
	
	neurons* nw=neurons_alloc(num_neurons,train_func_log,derivative_train_func_log);

	//number of points to fit
	int num_x=300;
	//Initial condition
	double x0=0;
	double y0=0.5;

	gsl_vector* vx=gsl_vector_alloc(num_x);

	for(int i=0;i<num_x;i++){
		double x=a+(b-a)*(i-1)/(num_x-1);
		gsl_vector_set(vx,i,x);
	}

	for(int i=0;i<nw->n;i++){
		gsl_vector_set(nw->data,0*nw->n+i,a+(b-a)*(i)/(nw->n-1));
		gsl_vector_set(nw->data,1*nw->n+i,1);
		gsl_vector_set(nw->data,2*nw->n+i,1);
	}

	neurons_train(nw,vx,Logistic_ODE,y0,x0);

	double dz=1.0/64;
	for(double z=a;z<=b;z+=dz){
		gsl_vector* FdF = gsl_vector_alloc(2);
		neurons_feed_forward(nw,z,FdF);
		printf("%g %g\n",z,gsl_vector_get(FdF,0));
		gsl_vector_free(FdF);

	}

	printf("\n\n");

////////////////////////////////////////////////////////////////////////////////
//Gaussian Differential equation:
////////////////////////////////////////////////////////////////////////////////

	//Number of hidden neurons to train
	num_neurons=10;
	
	neurons* nw_g=neurons_alloc(num_neurons,train_func_gauss,derivative_train_func_gauss);
	
	//Number of points to fit
	num_x=100;

	//Initial condition
	x0=0;
	y0=1;

	for(int i=0;i<nw_g->n;i++){
		gsl_vector_set(nw_g->data,0*nw_g->n+i,a+(b-a)*(i)/(nw_g->n-1));
		gsl_vector_set(nw_g->data,1*nw_g->n+i,1);
		gsl_vector_set(nw_g->data,2*nw_g->n+i,1);
	}

	neurons_train(nw_g,vx,Gaussian_ODE,y0,x0);

	for(double z=a;z<=b;z+=dz){
		gsl_vector* FdF_g = gsl_vector_alloc(2);
		neurons_feed_forward(nw_g,z,FdF_g);
		printf("%g %g\n",z,gsl_vector_get(FdF_g,0));
		gsl_vector_free(FdF_g);

	}
	
neurons_free(nw);
neurons_free(nw_g);
gsl_vector_free(vx);
return 0;
}