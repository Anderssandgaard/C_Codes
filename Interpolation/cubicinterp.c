#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

typedef struct {int n; double *x, *y, *b, *c, *d;} CubicSpline;

CubicSpline* CubicSpline_alloc(int n,double* x, double* y){

CubicSpline *s = (CubicSpline*)malloc(sizeof(CubicSpline));
 
s->b =(double*)malloc((n)*sizeof(double));
s->c =(double*)malloc((n-1)*sizeof(double));
s->d =(double*)malloc((n-1)*sizeof(double));
s->x =(double*)malloc((n)*sizeof(double));
s->y =(double*)malloc((n)*sizeof(double));
s->n = n;

for(int i=0; i<n; i++){
	s->x[i]=x[i];
	s->y[i]=y[i];
}

double p[n-1], h[n-1];

for(int i=0; i<n-1; i++){
	h[i]=x[i+1]-x[i];
	assert(h[i]>0);
	p[i]=(y[i+1]-y[i])/h[i];
}

double D[n], Q[n-1], B[n];
D[0]=2;
Q[0]=1;

for(int i=0; i<n-2; i++){
	D[i+1]=2*h[i]/h[i+1]+2;
	Q[i+1]=h[i]/h[i+1];
	B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]);
}

D[n-1]=2;
B[0]=3*p[0];
B[n-1]=3*p[n-2];

for(int i=0; i<n; i++){
	D[i]-=Q[i-1]/D[i-1];
	B[i]-=B[i-1]/D[i];
}

s->b[n-1]=B[n-1]/D[n-1];

for(int i=n-2;i>=0;i--){
	s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
}

for(int i=0; i<n; i++){
s->c[i]=(-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
s->d[i]=(s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
}


return s;
}


double CubicSpline_eval(CubicSpline *s, double Z){

assert(Z>=s->x[0] && Z<=s->x[s->n-1]);
int i=0, j=s->n-1;

while(j-i>1){
	int m=(i+j)/2;
	if(Z>s->x[m]){i=m;}
	else{j=m;}
	}

double h=Z-(s->x[i]);
return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));
}


double CubicSpline_eval_integ(CubicSpline *s, double Z){
assert(Z>=s->x[0] && Z<=s->x[s->n-1]);
int i=0, j=s->n-1;
while(j-i>1){
	int m=(i+j)/2;
	if(Z>s->x[m]){i=m;}
	else{j=m;}
	}
double Cubicintegral=0;

for(int k=1; k<=i-1; k++){
	double dxh=(s->x[k+1])-(s->x[k]);
	double I = dxh*s->y[k]+dxh*dxh*(0.5*s->b[k]+dxh*dxh*(1.0/3*s->c[k]+1.0/4*dxh*s->d[k]));
	Cubicintegral += I;
}
double h=Z-(s->x[i]);
Cubicintegral += h*s->y[i]+h*h*(0.5*s->b[i]+h*h*(1.0/3*s->c[i]+1.0/4*h*s->d[i]));
return Cubicintegral;
}

double CubicSpline_eval_deriv(CubicSpline *s, double Z){
assert(Z>=s->x[0] && Z<=s->x[s->n-1]);
int i=0, j=s->n-1;
while(j-i>1){
	int m=(i+j)/2;
	if(Z>s->x[m]){i=m;}
	else{j=m;}
	}
double h=Z-(s->x[i]);
double Cubicderivative =(s->b[i]+h*(2.0*s->c[i]+3.0*h*s->d[i]));
return Cubicderivative;
}



void CubicSpline_free(CubicSpline *s){
free(s->x);
free(s->y);
free(s->b);
free(s->c);
free(s->d);
free(s);
}
