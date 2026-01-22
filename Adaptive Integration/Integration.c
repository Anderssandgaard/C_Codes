#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include"matrixandvectors.h"
#include<assert.h>

double Adaptive_Integration24(double f(double x),double a,double b,double acc,double eps, double f2, double f3,int nrec,double *err){

assert(nrec<1e+6);

	double f1 = f(a+(b-a)/6);
	double f4 = f(a+5*(b-a)/6);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
	double q = (f1+f2+f3+f4)/4*(b-a);
	double tol = acc+eps*fabs(Q);
	*err = fabs(Q-q);

	if(*err<tol){
		return Q;
	}
	else{
		double Q1 = Adaptive_Integration24(f,a,(a+b)/2,acc/sqrt(2.),eps,f1,f2,nrec+1,err);
		double Q2 = Adaptive_Integration24(f,(a+b)/2,b,acc/sqrt(2.),eps,f3,f4,nrec+1,err);
		return Q1+Q2;
	}
}

double Adaptive_Integration(double f(double x),double a, double b, double acc, double eps, double *err){
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	int nrec = 0;
	return Adaptive_Integration24(f,a,b,acc,eps,f2,f3,nrec,err);
}

double Adaptive_Integration_Infinity(double f(double x),double a, double b, double acc, double eps,double *err){

int acheck=isinf(-a);
int bcheck=isinf(b);

	if(acheck && bcheck){
		double g(double t){
			return f(t/(1-t*t))*(1+t*t)/((1-t*t)*(1-t*t));
		};
		double a_new = -1, b_new = 1;
		double f2 = f(a_new+2*(b_new-a_new)/6);
		double f3 = f(a_new+4*(b_new-a_new)/6);
		int nrec = 0;
		return Adaptive_Integration24(g,a_new,b_new,acc,eps,f2,f3,nrec,err);
		
	}
	else if(acheck){
		double g(double t){
			return f(b-(1-t)/t)*1/(t*t);
		};
		double a_new = 0, b_new = 1;
		double f2 = f(a_new+2*(b_new-a_new)/6);
		double f3 = f(a_new+4*(b_new-a_new)/6);
		int nrec = 0;
		return Adaptive_Integration24(g,a_new,b_new,acc,eps,f2,f3,nrec,err);
		
	}

	else if(bcheck){
		double g(double t){
			return f(a+(1-t)/t)*1/(t*t);
		};
		double a_new = 0, b_new = 1;
		double f2 = f(a_new+2*(b_new-a_new)/6);
		double f3 = f(a_new+4*(b_new-a_new)/6);
		int nrec = 0;
		return Adaptive_Integration24(g,a_new,b_new,acc,eps,f2,f3,nrec,err);
		
	}

	else {
		
		double f2 = f(a+2*(b-a)/6);
		double f3 = f(a+4*(b-a)/6);
		int nrec = 0;
		return Adaptive_Integration24(f,a,b,acc,eps,f2,f3,nrec,err);
	}
}


double Clenshaw_Curtis(double f(double x),double a, double b, double acc, double eps, double *err){
double a_new=0;
double b_new=M_PI;

double Func_Clenshaw_Curtis(double theta){
	double x = 0.5*(a+b)+0.5*(a-b)*cos(theta);
	return f(x)*sin(theta)*(b-a)/2;
}
return Adaptive_Integration_Infinity(Func_Clenshaw_Curtis,a_new,b_new,acc,eps,err);
}


