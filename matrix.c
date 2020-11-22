#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

extern inline int m_index(const matrix *M, const int i, const int j) {
    return j*(M->rows) + i;
}

extern inline matrix* init_matrix(const int rows, const int cols) {
    matrix *M = (matrix *) calloc(1, sizeof(matrix));
    //printf("init %i %i\n", rows, cols);
    M->rows = rows;
    M->cols = cols;
    M->data = (double *) calloc(rows * cols, sizeof(double));
    return M;
}

extern inline matrix* copy_matrix(const matrix *M) {
    if (M == NULL) {
        return NULL;
    }

    matrix *M_copy = init_matrix(M->rows, M->cols);
    for (int j = 0; j < M->cols; j++) {
        for (int i = 0; i < M->rows; i++) {
            M_copy->data[m_index(M_copy, i, j)] = M->data[m_index(M, i, j)];
        }
    }
    // memcpy(M_copy->data, M->data, M->rows * M->cols * sizeof(double));

    return M_copy;
}

extern inline void free_matrix(matrix *M) {
    if (M != NULL) {
        free(M->data);
        free(M);
    }
}

matrix* mult(const matrix* restrict A, const matrix* restrict B) {
    if (A->cols != B->rows) {
        return NULL;
    }
    matrix* prod = init_matrix(A->rows, B->cols);
    for (int i = 0; i < prod->rows; i++) {
        for (int j = 0; j < prod->cols; j++) {
            double acc = 0.0;
            for (int k = 0; k < A->cols; k++) {
                
                acc += B->data[m_index(B, k, j)]*A->data[m_index(A, i, k)];
            }
            prod->data[m_index(prod, i, j)] = acc;
        }
    }

    return prod;
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

extern inline int sign(const double d) {
    return d > 0.0 ? 1 : (d < 0.0 ? -1 : 0);
}

void scalar_mult(const double c, matrix *M) {
    for (int j = 0; j < M->cols; j++) {
        for (int i = 0; i < M->rows; i++) {
            M->data[m_index(M, i, j)] *= c;
        }
    }
}

int add(matrix *A, const matrix *B) {
    if (A->rows != B->rows || A->cols != B-> cols) {
        return -1;
    }

    for (int j = 0; j < A->cols; j++) {
        for (int i = 0; i < A->rows; i++) {
            int index = m_index(A, i, j);
            A->data[index] += B->data[index];
        }
    }

    return 0;
}

int sub(matrix *A, const matrix *B) {
    if (A->rows != B->rows || A->cols != B-> cols) {
        return -1;
    }

    for (int j = 0; j < A->cols; j++) {
        for (int i = 0; i < A->rows; i++) {
            int index = m_index(A, i, j);
            A->data[index] -= B->data[index];
        }
    }

    return 0;
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

double norm(const matrix *v) {
    double n = 0.0;

    for (int i = 0; i < v->rows; i++) {
        double el = v->data[m_index(v, i, 0)];
        n += el*el;
    }

    return sqrt(n);
}

matrix* sub_matrix(const matrix *M, int i, int j, int k, int l) {
    if (M == NULL) {
        return NULL;
    }

    int rows = j - i + 1;
    int cols = l - k + 1;

    matrix *M_sub = init_matrix(rows, cols);

    for (int c = 0; c < cols; c++) {
        for (int r = 0; r < rows; r++) {
            M_sub->data[m_index(M_sub, r, c)] = M->data[m_index(M, i + r, k + c)];
        }
    }

    return M_sub;
}

void set_sub_matrix(matrix *M, int i, int j, int k, int l, const matrix *M_sub) {
    if (M == NULL || M_sub == NULL) {
        return;
    }

    int rows = j - i + 1;
    int cols = l - k + 1;

    for (int c = 0; c < cols; c++) {
        for (int r = 0; r < rows; r++) {
            M->data[m_index(M, i + r, k + c)] = M_sub->data[m_index(M_sub, r, c)];
        }
    }
}

QR* hqr(const matrix *M) {
    QR* qr = calloc(1, sizeof(QR));
    if (M == NULL) {
        qr->Q = NULL;
        qr->R = NULL;
        return qr;
    }

    int n = M->rows;
    int m = M->cols;

    matrix *Q = eye(n, n);
    matrix *R = copy_matrix(M);

    qr->Q = Q;
    qr->R = R;

    for (int j = 0; j < n; j++) {
        matrix *x = sub_matrix(R, j, n - 1, j, j);
        double normx = norm(x);
        double rjj = R->data[m_index(R, j, j)];
        int s = sign(rjj);

        double ul = rjj - s*normx;
        matrix *w = copy_matrix(x);
        scalar_mult(1.0/ul, w);
        w->data[m_index(w, 0, 0)] = 1.0;

        double tau = -s*ul/normx;

        matrix *R_sub = sub_matrix(R, j, n - 1, 0, m - 1);
        matrix *w_t = transpose(w);
        //printf("%i %i %i %i\n", w_t->rows, w_t->cols, R_sub->rows, R_sub->cols);   
        matrix *w_tTimesR_sub = mult(w_t, R_sub);
            
        matrix *R_diff = mult(w, w_tTimesR_sub);
        scalar_mult(tau, R_diff);
        sub(R_sub, R_diff);
        set_sub_matrix(R, j, n - 1, 0, m - 1, R_sub);

        matrix *Q_sub = sub_matrix(Q, 0, n - 1, j, m - 1);
        matrix *Q_subTimesw = mult(Q_sub, w);
        matrix *Q_diff = mult(Q_subTimesw, w_t);
        scalar_mult(tau, Q_diff);
        sub(Q_sub, Q_diff);
        set_sub_matrix(Q, 0, n - 1, j, m - 1, Q_sub);

        free_matrix(x);
        free_matrix(w);
        free_matrix(R_sub);
        free_matrix(w_t);
        free_matrix(w_tTimesR_sub);
        free_matrix(R_diff);
        free_matrix(Q_sub);
        free_matrix(Q_subTimesw);
        free_matrix(Q_diff);
    }

    return qr;
}

int main(int argc, char** argv) {
    setbuf(stdout, 0);
    matrix *B = init_matrix(3, 3);
    B->data[m_index(B, 0, 0)] = 1;
    B->data[m_index(B, 1, 0)] = 2;
    B->data[m_index(B, 2, 0)] = 3;

    B->data[m_index(B, 0, 1)] = 4;
    B->data[m_index(B, 1, 1)] = 5;
    B->data[m_index(B, 2, 1)] = 6;

    B->data[m_index(B, 0, 2)] = 7;
    B->data[m_index(B, 1, 2)] = 8;
    B->data[m_index(B, 2, 2)] = 9;

    QR *qr = hqr(B);
    print_matrix(B);
    matrix *Q_t = transpose(qr->Q);
    matrix *qrm = mult(qr->Q, Q_t);
    print_matrix(qrm);
    print_matrix(qr->R);

    free_matrix(B);
    free_matrix(qr->Q);
    free_matrix(qr->R);
    free(qr);
    free_matrix(qrm);
    free_matrix(Q_t);
   

   
    return 0;
}