#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include"matrixandvectors.h"


double fitfunc(int i, double arg){
   switch(i){
   case 0: return log(arg); break;
   case 1: return 1.0;   break;
   case 2: return arg;     break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}


int main(int argc, char** argv){
// Retrieve data
int n=atoi(argv[1]);
fprintf(stderr,"n=%i\n",n);
double e[n],s[n],k[n];
for(int i=0;i<n;i++){
	scanf("%lg %lg %lg",&(e[i]),&(s[i]),&(k[i]));
}

printf("---------------------------------------------------------------------\n");
printf("Least Squares Fit with SVD of data.txt to Function c1*log(x)+c2+c3*x:\n");
printf("---------------------------------------------------------------------\n");


int number_of_fit_function=3;


vector* x = vector_alloc(n);
vector* y = vector_alloc(n);
vector* dy = vector_alloc(n);
vector* c = vector_alloc(number_of_fit_function);
matrix* Sigma = matrix_alloc(number_of_fit_function,number_of_fit_function);
for(int i=0;i<n;i++){
	vector_set(x,i, e[i]);
	vector_set(y,i, s[i]);
	vector_set(dy,i,k[i]);
}

Sigma = Least_Square_fit_SVD(x,y,dy,number_of_fit_function,fitfunc,c);

vector* dc = vector_alloc(number_of_fit_function);
for(int i=0; i<number_of_fit_function;i++){
	vector_set(dc,i,sqrt(matrix_get(Sigma,i,i)));
}

double fitvalues(double x){
	double fx=0;
for(int j=0;j<number_of_fit_function;j++){
	fx+=vector_get(c,j)*fitfunc(j,x);
}
return fx;}


double fitvalues_minus_uncertainty(int j, double x){
	double fxmin=0;
	fxmin+=fitvalues(x)-fitfunc(j,x)*vector_get(dc,j);
return fxmin;}


double fitvalues_plus_uncertainty(int j, double x){
	double fxplus=0;
	fxplus+=fitvalues(x)+fitfunc(j,x)*vector_get(dc,j);
return fxplus;}

vector_print("Fitparameters c=\n",c);
matrix_print("\nCorrelation Matrix=",Sigma);
printf("\n\n");


double arg; 
double darg=(e[n-1]-e[0])/90;
for(int i=0;i<number_of_fit_function;i++){
arg=e[0]-darg;
do{
	printf("%g %g %g %g \n",arg,fitvalues(arg),fitvalues_plus_uncertainty(i,arg),fitvalues_minus_uncertainty(i,arg));
	arg+=darg;
}while(arg<e[n-1]+2);
printf("\n\n");
}

			
//FREE
vector_free(x);
vector_free(y);
vector_free(dy);
vector_free(c);
matrix_free(Sigma);
return 0;
}
