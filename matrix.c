#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/random.h>
#include <float.h>
#include <omp.h>
#include "matrix.h"

inline int m_index(const matrix *M, const int i, const int j) {
    return j*(M->rows) + i;
}

inline matrix* init_matrix(const int rows, const int cols) {
    matrix *M = (matrix *) calloc(1, sizeof(matrix));
    //printf("init %i %i\n", rows, cols);
    M->rows = rows;
    M->cols = cols;
    M->data = (double *) calloc(rows * cols, sizeof(double));
    return M;
}

inline matrix* copy_matrix(const matrix *M) {
    if (M == NULL) {
        return NULL;
    }

    matrix *M_copy = init_matrix(M->rows, M->cols);
    #pragma omp parallel for
    for (int j = 0; j < M->cols; j++) {
        for (int i = 0; i < M->rows; i++) {
            M_copy->data[m_index(M_copy, i, j)] = M->data[m_index(M, i, j)];
        }
    }
    // memcpy(M_copy->data, M->data, M->rows * M->cols * sizeof(double));

    return M_copy;
}

inline void free_matrix(matrix *M) {
    if (M != NULL) {
        free(M->data);
        free(M);
    }
}

matrix* mult(const matrix* restrict A, const matrix* restrict B) {
    if (A->cols != B->rows) {
        return NULL;
    }
    matrix* C = init_matrix(A->rows, B->cols);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->rows, B->cols, A->cols, 1.0, A->data, A->rows, B->data, B->rows, 0.0, C->data, A->rows);

    return C;
}

void print_matrix(const matrix *M) {
    if (M == NULL) {
        return;
    }
    printf("[ ");
    for (int i = 0; i < M->rows; i++) {
        if (i != 0) {
            printf("  ");
        }
        for (int j = 0; j < M->cols; j++) {
            printf("%lf ", M->data[m_index(M, i, j)]);
        }
        if (i != M->rows - 1) {
            printf("\n");
        }
    }
    printf("]\n\n");
}


matrix* eye(int rows, int cols) {
    matrix *M = init_matrix(rows, cols);
    for (int i = 0; i < rows && i < cols; i++) {
        M->data[m_index(M, i, i)] = 1.0;
    }

    return M;
}

matrix* transpose(const matrix *M) {
    matrix *M_t = init_matrix(M->cols, M->rows);
    for (int i = 0; i < M_t->rows; i++) {
        for (int j = 0; j < M_t->cols; j++) {
            M_t->data[m_index(M_t, i, j)] = M->data[m_index(M, j, i)];
        }
    }

    return M_t;
}

mask rand_mask(int rows, int cols) {
    unsigned char* m = calloc(rows*cols, sizeof(unsigned char));
    while (getrandom(m, rows*cols*sizeof(unsigned char), 0) < rows*cols*sizeof(bool)) { /* empty */ }
    for (int i = 0; i < rows*cols; i++) {
        m[i] = m[i] / 64;
    }

    return m;
}

void apply_mask(image_channels chns, mask m) {
    #pragma omp parallel for
    for (int j = 0; j < chns.R->cols; j++) {
        for (int i = 0; i < chns.R->rows; i++) {
            int index = m_index(chns.R, i, j);
            if (m[index]) {
                chns.R->data[index] = 255;
                chns.G->data[index] = 255;
                chns.B->data[index] = 255;
            }
        }
    }
}

image_channels image_to_channels(image img) {
    image_channels channels = {.R = init_matrix(img.height, img.width), .G = init_matrix(img.height, img.width), .B = init_matrix(img.height, img.width)};
    #pragma omp parallel for
    for (int i = 0; i < img.height; i++) {
        for (int j = 0; j < img.width; j++) {
            channels.R->data[m_index(channels.R, i, j)] = (double) img.rows[i][3*j];
            channels.G->data[m_index(channels.G, i, j)] = (double) img.rows[i][3*j + 1];
            channels.B->data[m_index(channels.B, i, j)] = (double) img.rows[i][3*j + 2];
        }
    }
    return channels;
}

image channels_to_image(image_channels channels) {
    png_bytepp row_pointers = calloc(channels.R->rows, sizeof(png_bytep));
    #pragma omp parallel for
    for (int i = 0; i < channels.R->rows; i++) {
        row_pointers[i] = calloc(3*channels.R->cols, sizeof(png_byte));
        for (int j = 0; j < channels.R->cols; j++) {
            row_pointers[i][3*j] = (png_byte) channels.R->data[m_index(channels.R, i, j)];
            row_pointers[i][3*j + 1] = (png_byte) channels.G->data[m_index(channels.G, i, j)];
            row_pointers[i][3*j + 2] = (png_byte) channels.B->data[m_index(channels.B, i, j)];
        }
    }

    image img = {.width = channels.R->cols, .height = channels.R->rows, .rows = row_pointers};
    return img;
}

matrix* ALS(const matrix* A, const int k) {
    matrix* W = eye(A->rows, k);
    matrix* Z = eye(A->cols, k);
    int jvpt[k];
    int rank;
    for (int i = 0; i < 30; i++) {
        matrix* A_cop = copy_matrix(A);
        matrix* W_cop = copy_matrix(W);
        // A = W, X = Z^T, B = A
        LAPACKE_dgelsy(LAPACK_COL_MAJOR, W_cop->rows, W_cop->cols, A->cols, W_cop->data, W_cop->rows, A_cop->data, A_cop->rows, jvpt, 1.e-8, &rank);
        free_matrix(W_cop);
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < A->cols; j++) {
                Z->data[m_index(Z, j, i)] = A_cop->data[m_index(A_cop, i, j)];
            }
        }
        free_matrix(A_cop);
        matrix* A_cop2 = transpose(A);
//print_matrix(A_cop2);
        matrix* Z_cop = copy_matrix(Z);
        // A = Z, X = W^T, B = A^T
        LAPACKE_dgelsy(LAPACK_COL_MAJOR, Z_cop->rows, Z_cop->cols, A_cop2->cols, Z_cop->data, Z->rows, A_cop2->data, A_cop2->rows, jvpt, 1.e-8, &rank);
        //print_matrix(A_cop2);
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < A->rows; j++) {
                W->data[m_index(W, j, i)] = A_cop2->data[m_index(A_cop2, i, j)];
            }
        }
        free_matrix(A_cop2);
        free_matrix(Z_cop);
    }
 //   print_matrix(W);
    //print_matrix(W);
    matrix* approx = init_matrix(A->rows, A->cols);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, W->rows, Z->rows, W->cols, 1.0, W->data, W->rows, Z->data, Z->rows, 0.0, approx->data, A->rows);
    free_matrix(W);
    free_matrix(Z);
    return approx;
}

matrix* ALS_masked(const matrix* A, mask m, const int k, const double beta, const int iterations) {
    matrix* W = eye(A->rows, k);
    matrix* Z = eye(A->cols, k);

    for (int iteration = 0; iteration < iterations; iteration++) {
        #pragma omp parallel for
        for (int i = 0; i < A->cols; i++) {
            int jvpt[k];
            int rank;
            int n = 0;
            for (int j = 0; j < A->rows; j++) {
                if (!m[m_index(A, j, i)]) {
                    n++;
                }
            }
            matrix* W_m = init_matrix(2*n, k);
            matrix* A_m = init_matrix(2*n, 1);
            n = 0;
            #pragma omp parallel for
            for (int j = 0; j < A->rows; j++) {
                if (!m[m_index(A, j, i)]) {
                    A_m->data[m_index(A_m, n, 0)] = A->data[m_index(A, j, i)];
                    for (int c = 0; c < W_m->cols; c++) {
                        W_m->data[m_index(W_m, n, c)] = W->data[m_index(W, j, c)];
                    }
                    n++;
                }
            }
            #pragma omp parallel for
            for (int r = n; r < 2*n; r++) {
                for (int j = 0; j < k; j++) {
                    if ((r-n) == j) {
                        W_m->data[m_index(W_m, r, j)] = beta;
                    }
                }
            }

            // A = W_m, X = Z_i, B = A_m
            int info;
            if ((info = LAPACKE_dgelsy(LAPACK_COL_MAJOR, W_m->rows, W_m->cols, A_m->cols, W_m->data, W_m->rows, A_m->data, A_m->rows, jvpt, 1.e-8, &rank)) != 0) {
                printf("%i\n", info);
                exit(info);
            }
            for (int j = 0; j < k; j++) {
                Z->data[m_index(Z, i, j)] = A_m->data[m_index(A_m, j, 0)];
            }
            free_matrix(W_m);
            free_matrix(A_m);
        }
        
        #pragma omp parallel for
        for (int i = 0; i < A->rows; i++) {
            int jvpt[k];
            int rank;
            int n = 0;
            for (int j = 0; j < A->cols; j++) {
                if (!m[m_index(A, i, j)]) {
                    n++;
                }
            }
            matrix* Z_m = init_matrix(2*n, k);
            matrix* A_m = init_matrix(2*n, 1);
            n = 0;
            #pragma omp parallel for
            for (int j = 0; j < A->cols; j++) {
                if (!m[m_index(A, i, j)]) {
                    A_m->data[m_index(A_m, n, 0)] = A->data[m_index(A, i, j)];
                    for (int c = 0; c < Z_m->cols; c++) {
                        Z_m->data[m_index(Z_m, n, c)] = Z->data[m_index(Z, j, c)];
                    }
                    n++;
                }
            }
            #pragma omp parallel for
            for (int r = n; r < 2*n; r++) {
                for (int j = 0; j < k; j++) {
                    if ((r - n) == j) {
                        Z_m->data[m_index(Z_m, r, j)] = beta;
                    }
                }
            }

            // A = Z_m, X = W_i, B = A_m
            if (LAPACKE_dgelsy(LAPACK_COL_MAJOR, Z_m->rows, Z_m->cols, A_m->cols, Z_m->data, Z_m->rows, A_m->data, A_m->rows, jvpt, 1.e-8, &rank) != 0) {
                exit(2);
            }
            for (int j = 0; j < k; j++) {
                W->data[m_index(W, i, j)] = A_m->data[m_index(A_m, j, 0)];
            }
            free_matrix(Z_m);
            free_matrix(A_m);
        }
    }
    matrix* approx = init_matrix(A->rows, A->cols);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, W->rows, Z->rows, W->cols, 1.0, W->data, W->rows, Z->data, Z->rows, 0.0, approx->data, A->rows);
    free_matrix(W);
    free_matrix(Z);
    return approx;
}

#define ARGS 6
int main(int argc, char** argv) {
    if (argc != ARGS + 1) {
        printf("Usage: %s SOURCE_IMG RANK BETA ITERATIONS THREADS OUTPUT_IMG\n", argv[0]);
        return 1;
    }
    image img = read_png(argv[1]);
    write_png(img, "orig.png");
    int rank = atoi(argv[2]);
    int iterations = atoi(argv[4]);
    int threads = atoi(argv[5]);
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
    double beta = strtod(argv[3], NULL);

    image_channels chns = image_to_channels(img);
    mask m = rand_mask(chns.R->rows, chns.R->cols);
    apply_mask(chns, m);
    image masked_img = channels_to_image(chns);
    write_png(masked_img, "masked.png");
    
    double start = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < 3; i++) {
        if (i == 0) {
            chns.R = ALS_masked(chns.R, m, rank, beta, iterations);
        } else if (i == 1) {
            chns.G = ALS_masked(chns.G, m, rank, beta, iterations);
        } else {
            chns.B = ALS_masked(chns.B, m, rank, beta, iterations);
        }
    }
    double end = omp_get_wtime();
    printf("Total time: %lfs\n", end - start);

    image new_img = channels_to_image(chns);
    

    write_png(new_img, argv[6]);
   
    return 0;
}   