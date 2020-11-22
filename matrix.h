#ifndef MATRIX_H
#define MATRIX_H

typedef struct matrix{
    int rows;
    int cols;

    // Stores data in column major order
    double* restrict data;
} matrix;

typedef struct QR {
    matrix* Q;
    matrix* R;
} QR;

// Return the array index for element starting at row i and col j (0-indexed)
inline int m_index(const matrix* M, int i, int j);

// Initializes a zero matrix with dimensions rows x cols
inline matrix* init_matrix(int rows, int cols);

// Returns a new matrix with the same values as M
inline matrix* copy_matrix(const matrix* M);

// Free the associated struct and data array for matrix M
inline void free_matrix(matrix* M);

// Multiplies matricies A (n_1 X k) and B (k X n_2). Returns a NULL pointer if dimensions don't match
matrix* mult(const matrix* A, const matrix* B);

void print_matrix(const matrix* A);

// Returns 1 if d is positive, 0 if d is zero, or -1 if d is negative
inline int sign(double d);

// Multiples a matrix M by scalar c in place
void scalar_mult(double c, matrix* M);

// Adds matrix B to A in place. Returns 0 if dimensions match, and -1 otherwise
int add(matrix *A, const matrix *B);

// Subtracts matrix B from A in place. Returns 0 if dimensions match, and -1 otherwise
int sub(matrix *A, const matrix *B);

// Returns a matrix with 1's on the main diagonal and 0 elsewhere
matrix* eye(int rows, int cols);

// Returns the transposed matrix M_t, where M_t(i, j) = M(j, i) 
matrix* transpose(const matrix *M);

// Returns the Euclidean norm (or magnitude) of the given column vector.
double norm(const matrix *v);

// Returns M(i:j, k:l) inclusive
matrix* sub_matrix(const matrix *M, int i, int j, int k, int l);

// Sets the values for M(i:j, k:l) to be the values of M_sub
void set_sub_matrix(matrix *M, int i, int j, int k, int l, const matrix *M_sub);

/* Returns the QR decomposition of M using householder transformations.
  That is, hqr returns matrices Q and R such that M = QR, Q is orthogonal and R is upper-triangular.
  If M has dimensions nxm, then Q is nxn and R is nxm.
*/
QR* hqr(const matrix* M);
#endif