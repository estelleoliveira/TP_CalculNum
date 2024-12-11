#include <stdio.h>
#include <cblas.h>
#include <lapacke.h>

int main() {
    printf("Test CBLAS and LAPACK\n");

    // Test BLAS: CBLAS example (ax = b)
    int n = 3;
    double alpha = 1.0;
    double A[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double x[3] = {1.0, 1.0, 1.0};
    double b[3];

    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, alpha, A, n, x, 1, 0.0, b, 1);

    printf("Result of CBLAS operation: ");
    for (int i = 0; i < n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n");

    // Test LAPACK: LAPACKE example (matrix inversion or other LAPACK functions)
    return 0;
}