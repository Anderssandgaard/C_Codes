#ifndef HAVE_SPLINE_H
#define HAVE_SPLINE_H

double linterp(int n, double *x, double *y, double Z);
double linterp_integ(int n, double *x, double *y, double Z);

typedef struct {int n; double *x, *y, *b, *c;} QuadSpline;
QuadSpline* QuadSpline_alloc(int n,double* x, double* y);
double QuadSpline_eval(QuadSpline *s, double Z);
double QuadSpline_eval_integ(QuadSpline *s, double Z);
double QuadSpline_eval_deriv(QuadSpline *s, double Z);
void QuadSpline_free(QuadSpline *s);

typedef struct {int n; double *x, *y, *b, *c, *d;} CubicSpline;
CubicSpline* CubicSpline_alloc(int n,double* x, double* y);
double CubicSpline_eval(CubicSpline *s, double Z);
double CubicSpline_eval_integ(CubicSpline *s, double Z);
double CubicSpline_eval_deriv(CubicSpline *s, double Z);
void CubicSpline_free(CubicSpline *s);
#endif
