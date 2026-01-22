#include<assert.h>
#include<stdio.h>
double linterp(int n, double *x, double *y, double Z){
assert(n>1 && Z>=x[0] && Z<=x[n-1]);
int i=0, j=n-1;
while(j-i>1){
	int m=(i+j)/2;
	if(Z>x[m]){i=m;}
	else{j=m;}
	}
return y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(Z-x[i]);
} 

double linterp_integ(int n, double *x, double *y, double Z){
assert(n>1 && Z>=x[0] && Z<=x[n-1]);
int i=0, j=n-1;
while(j-i>1){
	int m=(i+j)/2;
	if(Z>=x[m]){i=m;}
	else{j=m;}
	}

	double integral=0;

for(int k=1; k<i; k++){
	double dxk=x[k+1]-x[k];
	double bk=(y[k+1]-y[k])/dxk;
	double I = y[k]*dxk+0.5*bk*dxk*dxk;
	integral += I;
}

double dxZ = (Z-x[i]);
double bZ = (y[i+1]-y[i])/dxZ;
integral +=y[i]*dxZ+bZ*(0.5*dxZ*dxZ);
return integral;
}
