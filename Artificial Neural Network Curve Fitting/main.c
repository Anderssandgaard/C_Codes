#include<stdio.h>
#include<math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include"neurons.h"

double train_func(double x){return x*exp(-x*x);}

double func_to_interpolate(double x){return cos(5*x-1)*exp(-x*x);}

int main(){

	int num_neurons=20;
	
	neurons* nw=neurons_alloc(num_neurons,train_func);
	double a=-1,b=1;
	int num_x=20;
	gsl_vector* vx=gsl_vector_alloc(num_x);
	gsl_vector* vf=gsl_vector_alloc(num_x);
	for(int i=0;i<num_x;i++){
		double x=a+(b-a)*i/(num_x-1);
		double f=func_to_interpolate(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vf,i,f);
	}
	for(int i=0;i<nw->n;i++){
		gsl_vector_set(nw->data,0*nw->n+i,a+(b-a)*i/(nw->n-1));
		gsl_vector_set(nw->data,1*nw->n+i,1);
		gsl_vector_set(nw->data,2*nw->n+i,1);
	}

	neurons_train(nw,vx,vf);

	for(int i=0;i<vx->size;i++){
		double x=gsl_vector_get(vx,i);
		double f=gsl_vector_get(vf,i);
		printf("%g %g\n",x,f);
	}
	printf("\n\n");

	double dz=1.0/64;
	for(double z=a;z<=b;z+=dz){
		double y=neurons_feed_forward(nw,z);
		printf("%g %g\n",z,y);
	}
	
neurons_free(nw);
return 0;
}
