#include<stdio.h>
#include<math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include"neurons.h"

neurons* neurons_alloc(int n,double(*g)(double),double(*dg)(double)){
	neurons* nw = malloc(sizeof(neurons));
	nw->n=n;
	nw->g=g;
	nw->dg=dg;
	nw->data=gsl_vector_alloc(3*n);
	return nw;
}
void neurons_free(neurons* nw){
	gsl_vector_free(nw->data);
	free(nw);
}

void neurons_feed_forward(neurons* nw,double x,gsl_vector* FdF){
	double s, ds;
	for(int i=0;i<nw->n;i++){
		double a=gsl_vector_get(nw->data,0*nw->n+i);
		double b=gsl_vector_get(nw->data,1*nw->n+i);
		double w=gsl_vector_get(nw->data,2*nw->n+i);
		s+=(nw->g((x-a)/b)*w);
		ds+=(nw->dg((x-a)/b)*w/b);
	}
	gsl_vector_set(FdF,0,s);
	gsl_vector_set(FdF,1,ds);
}

void neurons_train(neurons* nw,gsl_vector* vx,double f(double x,double y),double y0,double x0){
	
	gsl_vector* FdF = gsl_vector_alloc(2);
	gsl_vector* FdF0 = gsl_vector_alloc(2);


	double delta(const gsl_vector* p, void* params){
	
		gsl_vector_memcpy(nw->data,p);
		double s=0;

		for(int i=0;i<vx->size;i++){
			double x=gsl_vector_get(vx,i);
			neurons_feed_forward(nw,x,FdF);
			double fun=f(x,gsl_vector_get(FdF,0));
			s+=pow(fabs(gsl_vector_get(FdF,1)-fun),2);
		}					
			neurons_feed_forward(nw,x0,FdF0);
		return s+vx->size*pow(fabs(gsl_vector_get(FdF0,0)-y0),2);
	}

	gsl_vector* p=gsl_vector_alloc(nw->data->size);
	gsl_vector_memcpy(p,nw->data);

	gsl_vector* step_size = gsl_vector_alloc(nw->data->size);
	gsl_vector_set_all(step_size, 0.1);

	gsl_multimin_function Func;
	Func.n = p->size;
	Func.f = delta;
	Func.params = NULL;

	gsl_multimin_fminimizer *s = 
	gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, Func.n);		
	gsl_multimin_fminimizer_set (s, &Func, p, step_size);
	
	int iter = 0, status;
	do{
		iter++;
		int iteration_status = gsl_multimin_fminimizer_iterate(s);
		if(iteration_status != 0)
		{
			fprintf(stderr, "No Improvement\n");
			break;
		}
		double acc = 0.0001;
		status = gsl_multimin_test_size(s->size, acc);
		if(status == GSL_SUCCESS) 
		{	
			fprintf(stderr, "-------------------------------------\n");
			fprintf(stderr, "Converged in %i iterations\n", iter);
			fprintf(stderr, "-------------------------------------\n");

		}

	}while( status == GSL_CONTINUE && iter < 1e+6);
	fprintf(stderr, "-------------------------------------\n");
	fprintf(stderr, "Neural Network Params\n");
	fprintf(stderr, "-------------------------------------\n");
	fprintf(stderr, "a b w\n");
	for (int i = 0; i < nw->n; i++)
	{
		fprintf(stderr, "%g \t %g \t %g \n", 
			gsl_vector_get(nw->data, 0*nw->n + i),
			gsl_vector_get(nw->data, 1*nw->n + i),
			gsl_vector_get(nw->data, 2*nw->n + i) );
	}

	gsl_vector_memcpy(nw->data, s->x); 
	gsl_vector_free(p);
	gsl_vector_free(FdF);
	gsl_vector_free(FdF0);
	gsl_vector_free(step_size);
	gsl_multimin_fminimizer_free(s);

}