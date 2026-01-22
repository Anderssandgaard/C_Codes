#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"matrixandvectors.h"
#include<assert.h>






void ODE_RK23_Step(double x, double h, vector* yx, void f(double x, vector* y, double* dydx),vector* yh, vector* error){

	int n=(*yx).size;
	int i;
	double k1[n],k2[n],k3[n],k4[n];
	vector* yt = vector_alloc(n);

	f(x,yx,k1);
	for(i=0;i<n;i++){
		vector_set(yt,i,vector_get(yx,i)+0.5*k1[1]*h);
	}
	f(x+0.5*h,yt,k2);
	for(i=0;i<n;i++){
		vector_set(yt,i,vector_get(yx,i)+0.75*k2[i]*h);
	}
	f(x+0.75*h,yt,k3);
	for(i=0;i<n;i++){
		vector_set(yh,i,vector_get(yx,i)+(2./9*k1[i]+1./3*k1[i]+4./9*k1[i])*h);
	}
	f(x+h,yh,k4);
	for(i=0;i<n;i++){
		vector_set(yt,i,vector_get(yx,i)+(7./24*k1[i]+1./4*k2[i]+1./3*k3[i]+1./8*k4[i])*h);
		vector_set(error,i,vector_get(yh,i)-vector_get(yt,i));
	}
vector_free(yt);
}



int ODE_Driver(void f(double x, vector* yx, double* dydx), matrix* ylist, vector* xlist, double b, double h, double acc, double eps, int step_max){
	
    int step=0;
    int n=(*ylist).size2;
	double x; 
    double err;
    double norm_y;
    double tol;
    double a=vector_get(xlist, 0);
	vector* yh=vector_alloc(n);
	vector* dy=vector_alloc(n);
	vector* y=vector_alloc(n);


	while(vector_get(xlist,step)<b){

		x=vector_get(xlist, step);
		y=matrix_get_row(ylist,step);

		if(x+h>b){
			h=b-x;
		}		

		ODE_RK23_Step(x,h,y,f,yh,dy);

		err=sqrt(vector_dot_product(dy,dy));
		norm_y=sqrt(vector_dot_product(yh,yh));
		tol=(norm_y*eps+acc)*sqrt(h*1.0/(b-a));

		if(err<tol){
			step++;
			if(step>step_max-1){
				fprintf(stderr, "\n\nMax num of steps reached.\n\n");
				break;
			}
			vector_set(xlist, step, x+h);
			for(int i=0;i<n;i++){
			matrix_set(ylist, step,i,vector_get(yh,i));
			}
		}
		if(err>0){
			h*=pow(tol/err,0.25)*0.95;
		}
		else{
			h*=2;
		}
	}
	vector_free(yh);
	vector_free(dy);
	return step+1;
}


double ODE_Integral(void f(double x, vector* yx, double* dydx),double a,vector* init, double b,int numfunc, double h, double acc, double eps, int step_max){
 
matrix* ylist=matrix_alloc(step_max, numfunc);
vector* xlist=vector_alloc(step_max);

vector_set(xlist, 0, a);

for(int i=0;i<numfunc;i++){
matrix_set(ylist,0,i,vector_get(init,i));
}

ODE_Driver(f,ylist,xlist,b,h,acc,eps,step_max);
double s=0;
for(int i=1;i<(*ylist).size1;i++){
    s+=matrix_get(ylist,i,0)*(vector_get(xlist,i)-vector_get(xlist,i-1));
}




matrix_free(ylist);
vector_free(xlist);
return s;
}















