#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include"matrixandvectors.h"


void System_of_Eqn(vector* p, vector* fx){
	int A = 10000;
	double x=vector_get(p,0);
	double y=vector_get(p,1);
	vector_set(fx,0, A*x*y-1);
	vector_set(fx,1, exp(-x)+exp(-y)-1.-1.0/A);
	}

void System_of_Eqn_Jacobian(vector* p, matrix* J){
	int A = 10000;
	double x=vector_get(p,0);
	double y=vector_get(p,1);
	matrix_set(J,0,0,A*y);
	matrix_set(J,0,1,A*x);
	matrix_set(J,1,0,-exp(-x));
	matrix_set(J,1,1,-exp(-y));
	}

void Rosenbrock(vector* p,vector* fx){
	double x=vector_get(p,0);
	double y=vector_get(p,1);
	vector_set(fx,0, -2*(1-x)+400*(y-x*x)*(-1)*x);
	vector_set(fx,1, 200*(y-x*x));
	}

void Rosenbrock_Jacobian(vector* p, matrix* J){
	double x=vector_get(p,0);
	double y=vector_get(p,1);
	matrix_set(J,0,0,2-400*(y-3*x*x));
	matrix_set(J,0,1,-400*x);
	matrix_set(J,1,0,-400*x);
	matrix_set(J,1,1,200);
	}

void Himmelblau(vector* p,vector* fx){
	double x=vector_get(p,0);
	double y=vector_get(p,1);
	vector_set(fx,0, 4*x*(x*x+y-11)+2*(x+y*y-7));
	vector_set(fx,1, 2*(x*x+y-11)+4*y*(x+y*y-7));
	}

void Himmelblau_Jacobian(vector* p, matrix* J){
	double x=vector_get(p,0);
	double y=vector_get(p,1);
	matrix_set(J,0,0,4*(3*x*x+y-11)+2);
	matrix_set(J,0,1,4*x+4*y);
	matrix_set(J,1,0,4*x+4*y);
	matrix_set(J,1,1,2+4*(x+3*y*y-7));
	}

int main(){

vector* x=vector_alloc(2);
vector* fx=vector_alloc(2);
vector* numcalls=vector_alloc(2);

printf("---------------------------------------------------------------------\n");
printf("ROOT FINDING FOR VARIOUS FUNCTIONS \n(SEE OUT.GSL.TXT FOR THE GSL GENERATED VALUES:\n");
printf("GSL w.o. Jacobian for System of Eqn: 165 Steps\n");
printf("GSL w.o. Jacobian for Rosenbrock' Function: 132 Steps\n");
printf("GSL w.o. Jacobian for Himmelblau's Function: 7 Steps\n");
printf("---------------------------------------------------------------------\n\n");


//SYSTEM OF EQUATIONS (Normal line search vs. w. Jacobian vs. w. Jacobain and Refined polynomial linesearch)


printf("---------------------------------------------------------------------\n");
printf("Roots for System of Eqn A*x*y=1, exp(-x)+exp(-y)=1+1/A, A=10000:\n");
printf("---------------------------------------------------------------------\n");

		vector_set(x,0,0);
		vector_set(x,1,1);

		vector_print("\nInitial Vector x (Starting Point for line search)\n:",x);
		System_of_Eqn(x,fx);
		vector_print("\nf_(System_of_Eqn)(x):\n",fx);


//W. Normal Linesearch:
		
numcalls = Find_Roots_Newton(System_of_Eqn,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		System_of_Eqn(x,fx);
		vector_print("\nf(x_root):\n",fx);

// W. Jacobian

		vector_set(x,0,0);
		vector_set(x,1,1);

numcalls = Find_Roots_Newton_Jacobian(System_of_Eqn,System_of_Eqn_Jacobian,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch with Analytic Jacobian:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		System_of_Eqn(x,fx);
		vector_print("\nf(x_root):\n",fx);

// W. Jacobian and Refined Line Search
		
		vector_set(x,0,0);
		vector_set(x,1,1);

numcalls = Find_Roots_Newton_Jacobian_Refined(System_of_Eqn,System_of_Eqn_Jacobian,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch with Analytic Jacobian and Refined Line Search:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		System_of_Eqn(x,fx);
		vector_print("\nf(x_root):\n",fx);



//ROSENBROCK (Normal line search vs. w. Jacobian vs. w. Jacobain and Refined polynomial linesearch)
printf("\n\n---------------------------------------------------------------------\n");
printf("---------------------------------------------------------------------\n");
printf("Roots for Rosenbrock's Function (1-x)^2+100(y-x^2)^2:\n");
printf("---------------------------------------------------------------------\n");
printf("---------------------------------------------------------------------\n");

		vector_set(x,0,-2);
		vector_set(x,1,-2);

		vector_print("\nInitial Vector x (Starting Point for line search):\n",x);
		Rosenbrock(x,fx);
		vector_print("\nf_(Rosenbrock)(x):\n",fx);


//W. Normal Linesearch:

	

numcalls = Find_Roots_Newton(Rosenbrock,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		Rosenbrock(x,fx);
		vector_print("\nf(x_root):\n",fx);

// W. Jacobian

		vector_set(x,0,-2);
		vector_set(x,1,-2);

numcalls = Find_Roots_Newton_Jacobian(Rosenbrock,Rosenbrock_Jacobian,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch with Analytic Jacobian:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		Rosenbrock(x,fx);
		vector_print("\nf(x_root):\n",fx);

// W. Jacobian and Refined Line Search
		
		vector_set(x,0,-2);
		vector_set(x,1,-2);

numcalls = Find_Roots_Newton_Jacobian_Refined(Rosenbrock,Rosenbrock_Jacobian,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch with Analytic Jacobian and Refined Line Search:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		Rosenbrock(x,fx);
		vector_print("\nf(x_root):\n",fx);


//HIMMELBLAU (Normal line search vs. w. Jacobian vs. w. Jacobain and Refined polynomial linesearch)

printf("\n\n---------------------------------------------------------------------\n");
printf("---------------------------------------------------------------------\n");
printf("Roots for Himmelblau's Function (x^2+y-11)^2+(x+y^2-7)^2:\n");
printf("---------------------------------------------------------------------\n");
printf("---------------------------------------------------------------------\n");


		vector_set(x,0,0);
		vector_set(x,1,0);

		vector_print("\nInitial Vector x (Starting Point for line search):\n",x);
		Himmelblau(x,fx);
		vector_print("\nf_(Himmelblau)(x):\n",fx);


//W. Normal Linesearch:
		

numcalls = Find_Roots_Newton(Himmelblau,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		Himmelblau(x,fx);
		vector_print("\nf(x_root):\n",fx);

// W. Jacobian

		vector_set(x,0,0);
		vector_set(x,1,0);

numcalls = Find_Roots_Newton_Jacobian(Himmelblau,Himmelblau_Jacobian,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch with Analytic Jacobian:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		Himmelblau(x,fx);
		vector_print("\nf(x_root):\n",fx);

// W. Jacobian and Refined Line Search

		vector_set(x,0,0);
		vector_set(x,1,0);

numcalls = Find_Roots_Newton_Jacobian_Refined(Himmelblau,Himmelblau_Jacobian,x,1e-6,1e-3);

printf("---------------------------------------------------------------------\n");
printf("For Normal Linesearch with Analytic Jacobian and Refined Line Search:\n");
printf("---------------------------------------------------------------------\n");

		vector_print("\n1. Column = Number of steps. 2.Column = number of function calls:\n",numcalls);
		vector_print("\nsolution x_root:\n",x);
		Himmelblau(x,fx);
		vector_print("\nf(x_root):\n",fx);



//FREE ALLOCTIONS
vector_free(x);
vector_free(fx);
vector_free(numcalls);

return 0;
}
