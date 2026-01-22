#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include"neurons.h"

neurons* neurons_alloc(int n,double(*f)(double, double)){
	neurons* nw = malloc(sizeof(neurons));
	nw->n=n;
	nw->f=f;
	nw->data=gsl_vector_alloc(5*n);
	return nw;
}
void neurons_free(neurons* nw){
	gsl_vector_free(nw->data);
	free(nw);
}

double neurons_feed_forward(neurons* nw,double x, double y){
	double s=0;
	for(int i=0;i<nw->n;i++){
		double a=gsl_vector_get(nw->data,0*nw->n+i);
		double b=gsl_vector_get(nw->data,1*nw->n+i);	
		double c=gsl_vector_get(nw->data,2*nw->n+i);
		double d=gsl_vector_get(nw->data,3*nw->n+i);
		double w=gsl_vector_get(nw->data,4*nw->n+i);
		s+=nw->f((x+a)/b,(y+c)/d)*w;
	}
	return s;
}

void neurons_train(neurons* nw,gsl_vector* vx,gsl_vector* vy,gsl_vector* vf){
	double delta(const gsl_vector* p, void* params){
		gsl_vector_memcpy(nw->data,p);
		double s=0;
		for(int i=0;i<vx->size;i++){
			double x=gsl_vector_get(vx,i);
			double y=gsl_vector_get(vy,i);
			double f=gsl_vector_get(vf,i);
			double yy=neurons_feed_forward(nw,x,y);
			s+=fabs(yy-f);
		}
		return s/vx->size;
	}
	gsl_vector* p=gsl_vector_alloc(nw->data->size);
	gsl_vector_memcpy(p,nw->data);

	gsl_vector* step_size = gsl_vector_alloc(nw->data->size);
	gsl_vector_set_all(step_size, 0.1);

	gsl_multimin_function F;
	F.n = p->size;
	F.f = delta;
	F.params = NULL;

	gsl_multimin_fminimizer *s = 
	gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, F.n);		
	gsl_multimin_fminimizer_set (s, &F, p, step_size);
	
	int iter = 0, status;
	do{
		iter++;
		int iteration_status = gsl_multimin_fminimizer_iterate(s);
			if(iteration_status != 0)
			{
			fprintf(stderr, "No Improvement\n");
			break;
			}

		double acc = 0.001;
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
	fprintf(stderr, "a b c d w\n");
	
			for (int i = 0; i < nw->n; i++)
			{
			fprintf(stderr, "%g \t %g \t %g \t %g \t %g \n", 
			gsl_vector_get(nw->data, 0*nw->n + i),
			gsl_vector_get(nw->data, 1*nw->n + i),
			gsl_vector_get(nw->data, 2*nw->n + i),
			gsl_vector_get(nw->data, 3*nw->n + i),
			gsl_vector_get(nw->data, 4*nw->n + i) );
			}

	gsl_vector_memcpy(nw->data, s->x); 
	gsl_vector_free(p);
	gsl_vector_free(step_size);
	gsl_multimin_fminimizer_free(s);

}
