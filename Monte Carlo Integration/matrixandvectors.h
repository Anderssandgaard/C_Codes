#ifndef HAVE_MATRIXANDVECTORS_H
#define HAVE_MATRIXANDVECTORS_H

//VECTORS//
typedef struct { int size; double *data;} vector;

vector* vector_alloc(int n);

void vector_set(vector* v, int i, double vi);

double vector_get(const vector* v, int i);

void vector_print(const char* s, vector* v);

double  vector_dot_product (vector* u, vector* v);

void vector_add (vector* v, vector* u, double alpha);

void vector_set_unit(vector* v,int j);

vector* vector_copy(const vector* v);

void vector_free(vector* v);

void vector_scale(vector* v, double r);


//MATRICES//
typedef struct { int size1, size2; double *data;} matrix;

matrix* matrix_alloc(int n, int m);

void matrix_set(matrix* A, int i, int j, double aij);

double matrix_get(matrix* A, int i, int j);

void matrix_print(const char* s, matrix* A);

vector* matrix_get_column(matrix* A,int i);

vector* matrix_get_row(matrix* A,int i);

void matrix_set_column(matrix* A,int i, vector* v);

void matrix_set_row(matrix* A,int i, vector* v);

matrix* matrix_transpose(matrix* A);

matrix* matrix_multiplication(matrix* A,char TA,matrix* B,char TB);

vector* matrix_vector_multiplication(matrix* A,const char T,const vector* x);

vector* vector_matrix_multiplication(const vector* x,matrix* A,const char T);

matrix* matrix_copy(matrix* A);

void matrix_free(matrix* A);

void matrix_identity(matrix* A);

matrix* vector_outer_product(vector* x,vector* y);

void matrix_add(matrix* A, matrix* B);

void matrix_scale(matrix* A,double a);

// QR DECOMPOSITION: GRAM SCHMIDT AND GIVENS METHOD.
void Gram_Schmidt_QR_Decomposition(matrix* A, matrix* R);

void BackSubstitution(matrix* R,vector* x);

void Gram_Schmidt_QR_Solve(matrix* Q,matrix* R,vector* b, vector* x);

void Gram_Schmidt_QR_Inverse(matrix* Q,matrix* R,matrix* B);

void Givens_QR_Decomposition(matrix* A);

void Givens_QR_Solve(matrix* QR, vector* b);

void Givens_QR_QTb(matrix* QR, vector* b);

void Givens_QR_Inverse(matrix* QR,matrix* B);

// MATRIX DIAGONALISATION

int Jacobi_Diagonalisation(matrix* A ,vector* x, matrix* V);

int Jacobi_Diagonalisation_Classical(matrix* A ,vector* x, matrix* V);

int Jacobi_Diagonalisation_eigbyeig(matrix* A ,vector* x, matrix* V,int order,int num_of_eig);

void Singular_Value_Decomposition(matrix*A,matrix* S,matrix* V);


//LEAST SQUARES METHODS

matrix* Least_Square_fit(vector* x, vector* y, vector* dy, int numfunc, double func(int n, double x), vector* c);

matrix* Least_Square_fit_SVD(vector* x, vector* y, vector* dy, int numfunc, double func(int n, double x), vector*c);




// ROOT FINDING

vector* Find_Roots_Newton(void func(vector* x, vector* fx),vector* x,double dx,double eps);

vector* Find_Roots_Newton_Jacobian(void func(vector* x, vector* fx),void Jacobian(vector* x, matrix* J),vector* x,double dx,double eps);

vector* Find_Roots_Newton_Jacobian_Refined(void func(vector* x, vector* fx),void Jacobian(vector* x, matrix* J),vector* x,double dx,double eps);


// MINIMISATION

vector* Newton_Minimisation(double func(vector* x, vector* Gradient, matrix* Hessian),vector* x,double eps);

vector* Newton_Minimisation_Broyden(double f(vector* x, vector* Gradient), vector* x, double eps);

double** make_simplex(double* x, int d);

void print_simplex(double** simplex,int d);

double** alloc_matrix(int rows, int cols);

void free_matrix(double** matrix,int rows, int cols);

int downhill_simplex(double F(double*),double** simplex, int d, double simplex_size_goal);



// ODE

void ODE_RK23_Step(double x, double h, vector* y, void f(double x, vector* y, double* dydx),vector* yh, vector* error);

int ODE_Driver(void f(double x, vector* yx, double* dydx), matrix* ylist, vector* xlist, double b, double h, double acc, double eps, int step_max);

double ODE_Integral(void f(double x, vector* yx, double* dydx),double a,vector* init, double b,int numfunc, double h, double acc, double eps, int step_max);



// ADAPTIVE INTEGRATION

double Adaptive_Integration(double f(double),double a, double b, double acc, double eps, double *err);

double Adaptive_Integration24(double f(double x),double a,double b,double acc,double eps, double f2, double f3,int nrec,double *err);

double Adaptive_Integration_Infinity(double f(double x),double a, double b, double acc, double eps,double *err);


double Clenshaw_Curtis(double f(double x),double a, double b, double acc, double eps, double *err);


// MONTE CARLO INTEGRATION

void Monte_Carlo_Plain(double f(vector* x), vector* a, vector* b, int dim, int N, double *res, double *err);		

double Adaptive_2D_Integrator(double f(double x, double y),double c(double x), double d(double x),double a, double b,double acc, double eps, double *err);










#endif	

