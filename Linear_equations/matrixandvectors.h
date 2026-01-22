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

void matrix_set_column(matrix* A,int i, vector* v);

matrix* matrix_transpose(matrix* A);

matrix* matrix_multiplication(matrix* A,char TA,matrix* B,char TB);

vector* matrix_vector_multiplication(matrix* A,const char T,const vector* x);

matrix* matrix_copy(matrix* A);

void matrix_free(matrix* A);

// QR DECOMPOSITION: GRAM SCHMIDT AND GIVENS METHOD.
void Gram_Schmidt_QR_Decomposition(matrix* A, matrix* R);

void BackSubstitution(matrix* R,vector* x);

void Gram_Schmidt_QR_Solve(matrix* Q,matrix* R,vector* b, vector* x);

void Gram_Schmidt_QR_Inverse(matrix* Q,matrix* R,matrix* B);

void Givens_QR_Decomposition(matrix* A);

void Givens_QR_Solve(matrix* QR, vector* b);

void Givens_QR_QTb(matrix* QR, vector* b);

void Givens_QR_Inverse(matrix* QR,matrix* B);

#endif	

