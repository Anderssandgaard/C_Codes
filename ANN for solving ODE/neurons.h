#include<gsl/gsl_vector.h>
#ifndef HAVE_NEURONS
#define HAVE_NEURONS
typedef struct {
	int n;
	double(*g)(double);
	double(*dg)(double);
	gsl_vector* data;
	} neurons;
neurons* neurons_alloc(int n,double(*g)(double),double(*dg)(double));
void neurons_free(neurons* nw);
void neurons_feed_forward(neurons* nw,double x,gsl_vector* FdF);
void neurons_train(neurons* nw,gsl_vector* xlist,double(*f)(double,double),double y0,double x0);
#endif