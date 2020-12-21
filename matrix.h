#ifndef MATRIX_H
#define MATRIX_H
#include <lapacke.h>
#include <cblas.h>
#include <stdbool.h>
#include "png.h"

typedef struct matrix{
    int rows;
    int cols;

    // Stores data in column major order
    double* restrict data;
} matrix;

typedef bool* mask;

typedef struct {
    matrix* R;
    matrix* G;
    matrix* B;
} image_channels;

// Return the array index for element starting at row i and col j (0-indexed)
extern inline int m_index(const matrix* M, const int i, const int j);

// Initializes a zero matrix with dimensions rows x cols
extern inline matrix* init_matrix(const int rows, const int cols);

// Returns a new matrix with the same values as M
extern inline matrix* copy_matrix(const matrix* M);

// Free the associated struct and data array for matrix M
extern inline void free_matrix(matrix* M);

// Prints each row of matrix A sequentially
void print_matrix(const matrix* A);

// Returns a matrix with 1's on the main diagonal and 0 elsewhere
matrix* eye(int rows, int cols);

// Returns the transposed matrix M_t, where M_t(i, j) = M(j, i) 
matrix* transpose(const matrix *M);

// Randomly generated a mask, where m[(i, j)] == true means that the pixel at (i, j) is unknown/noise
mask rand_mask(int rows, int cols);

// Set all pixels (i, j) where m[(i, j)] == true to white
void apply_mask(image_channels chns, mask m);

// Combine three matrices corrresponding to the RGB channels to an image
image channels_to_image(image_channels channels);

// Split an image into three matrices corresponding to the RGB channels
image_channels image_to_channels(image img);

// Generate a rank k approximation of matrix A using alternating least squares
matrix* ALS(const matrix* A, const int k);

// Generate a rank k approximation of matrix A with missing values described by mask.
matrix* ALS_masked(const matrix* A, mask m, const int k, const double beta, const int iterations);
#endif