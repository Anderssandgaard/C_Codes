#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>



int main(int argc, char** argv){


int n=atoi(argv[1]);
fprintf(stderr,"n=%i\n",n);
double e[n],s[n];
for(int i=0;i<n;i++){
	scanf("%lg %lg",&(e[i]),&(s[i]));
}


double dx=0.1;
{ double xi, yi	;
    gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_linear, n);

    gsl_spline_init (spline, e, s, n);

    for (xi = e[0]; xi < e[9]; xi += dx)
      {
        yi = gsl_spline_eval (spline, xi, acc);
        printf ("%g %g\n", xi, yi);
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }




return 0;
}