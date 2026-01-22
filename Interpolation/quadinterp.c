#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

typedef struct {int n; double *x, *y, *b, *c;} QuadSpline;

QuadSpline* QuadSpline_alloc(int n,double* x, double* y){
QuadSpline *s = (QuadSpline*)malloc(sizeof(QuadSpline));
 
s->b =(double*)malloc((n-1)*sizeof(double));
s->c =(double*)malloc((n-1)*sizeof(double));
s->x =(double*)malloc((n)*sizeof(double));
s->y =(double*)malloc((n)*sizeof(double));
s->n = n;

for(int i=0; i<n; i++){
	s->x[i]=x[i];
	s->y[i]=y[i];
}

int i;
double p[n-1], h[n-1];
for(i=0; i<n-1; i++){
	h[i]=x[i+1]-x[i];
	assert (h[i]>0);
	p[i]=(y[i+1]-y[i])/h[i];
}
s->c[0]=0;
for(i=0; i<n-2; i++){
	s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
}
s->c[n-2]/=2;
for(i=n-3; i>=0; i--){
	s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
}
for(i=0; i<n-1; i++){
	s->b[i]=p[i]-s->c[i]*h[i];
}

return s;
}


double QuadSpline_eval(QuadSpline *s, double Z){
assert(Z>=s->x[0] && Z<=s->x[s->n-1]);
int i=0, j=s->n-1;
while(j-i>1){
	int m=(i+j)/2;
	if(Z>s->x[m]){i=m;}
	else{j=m;}
	}
double h=Z-(s->x[i]);
return s->y[i]+h*(s->b[i]+h*s->c[i]);
}


double QuadSpline_eval_integ(QuadSpline *s, double Z){
assert(Z>=s->x[0] && Z<=s->x[s->n-1]);
int i=0, j=s->n-1;
while(j-i>1){
	int m=(i+j)/2;
	if(Z>s->x[m]){i=m;}
	else{j=m;}
	}
double Quadintegral=0;

for(int k=1; k<=i-1; k++){
	double dxh=(s->x[k+1])-(s->x[k]);
	double I = dxh*(s->y[k])+0.5*dxh*dxh*(s->b[k])+1.0/3*(dxh*dxh*dxh)*(s->c[k]);
	Quadintegral += I;
}
double h=Z-(s->x[i]);
Quadintegral += h*(s->y[i])+0.5*h*h*(s->b[i])+1.0/3*(h*h*h)*(s->c[i]);
return Quadintegral;
}

double QuadSpline_eval_deriv(QuadSpline *s, double Z){
assert(Z>=s->x[0] && Z<=s->x[s->n-1]);
int i=0, j=s->n-1;
while(j-i>1){
	int m=(i+j)/2;
	if(Z>s->x[m]){i=m;}
	else{j=m;}
	}
double h=Z-(s->x[i]);
double Quadderivative =(s->b[i]+2.0*h*(s->c[i]));
return Quadderivative;
}



void QuadSpline_free(QuadSpline *s){
free(s->x);
free(s->y);
free(s->b);
free(s->c);
free(s);
}
