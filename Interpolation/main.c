#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include"spline.h"


int main(int argc, char** argv){

// Retrieve data
int n=atoi(argv[1]);
fprintf(stderr,"n=%i\n",n);
double e[n],s[n];
for(int i=0;i<n;i++){
	scanf("%lg %lg",&(e[i]),&(s[i]));
}
//

//Allocation
QuadSpline* Quad = QuadSpline_alloc(n,e,s);
CubicSpline* Cubic = CubicSpline_alloc(n,e,s);
//


// GSL Engines
gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();

gsl_spline *spline_lin
      = gsl_spline_alloc (gsl_interp_linear, n);
gsl_spline_init (spline_lin, e, s, n);

gsl_spline *spline_cubic
      = gsl_spline_alloc (gsl_interp_cspline, n);
gsl_spline_init (spline_cubic, e, s, n);
//

//Get data and print
double dx=0.01;
for(double x=e[0]; x<e[n-1]; x+=dx){

	double spline_res = linterp(n,e,s,x);
	double spline_res_integ = linterp_integ(n,e,s,x);
	
	double spline_quad_res = QuadSpline_eval(Quad,x);
	double spline_quad_integ_res = QuadSpline_eval_integ(Quad,x);
	double spline_quad_deriv_res = QuadSpline_eval_deriv(Quad,x);
	
	double spline_cubic_res = CubicSpline_eval(Cubic,x);
	double spline_cubic_integ_res = CubicSpline_eval_integ(Cubic,x);
	double spline_cubic_deriv_res = CubicSpline_eval_deriv(Cubic,x);
	
	double gsl_spline_lin_res = gsl_spline_eval (spline_lin, x, acc);
    double gsl_spline_lin_integ_res = gsl_spline_eval_integ (spline_lin, 0, x, acc);	
	
	double gsl_spline_cubic_res = gsl_spline_eval (spline_cubic, x, acc);
    double gsl_spline_cubic_integ_res = gsl_spline_eval_integ (spline_cubic, 0, x, acc);
	double gsl_spline_cubic_deriv_res = gsl_spline_eval_deriv (spline_cubic, x, acc);

		printf("%g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
							x,
							spline_res,
							spline_res_integ,
								spline_quad_res,
								spline_quad_integ_res,
								spline_quad_deriv_res,
									gsl_spline_lin_res,
									gsl_spline_lin_integ_res,
										gsl_spline_cubic_res,
										gsl_spline_cubic_integ_res,
										gsl_spline_cubic_deriv_res,
											spline_cubic_res,
											spline_cubic_integ_res,
											spline_cubic_deriv_res);
				}

// Clean
QuadSpline_free(Quad);
gsl_spline_free (spline_lin);
gsl_spline_free (spline_cubic);
gsl_interp_accel_free (acc);
return 0;
}
