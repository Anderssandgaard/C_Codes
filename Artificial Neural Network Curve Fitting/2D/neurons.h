#include<gsl/gsl_vector.h>
#ifndef HAVE_NEURONS
#define HAVE_NEURONS
typedef struct {
	int n;
	double(*f)(double,double);
	gsl_vector* data;
	} neurons;
neurons* neurons_alloc(int n,double(*f)(double,double));
void neurons_free(neurons* nw);
double neurons_feed_forward(neurons* nw,double x,double y);
void neurons_train(neurons* nw,gsl_vector* xlist,gsl_vector* ylist,gsl_vector* flist);
#endif