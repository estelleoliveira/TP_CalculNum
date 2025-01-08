/**********************************************/
/* lib_poisson1D.h                            */
/* Header for Numerical library developed to  */ 
/* solve 1D Poisson problem (Heat equation)   */
/**********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "atlas_headers.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int *la, int *kv);
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int *la, int *kv);
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1);
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1);
void set_grid_points_1D(double* x, int* la);
double relative_forward_error(double* x, double* y, int* la);
void write_GB2AIJ_operator_poisson1D(double* AB, int *la, char* filename);
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int *la, char* filename);
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename);
void write_vec(double* vec, int* la, char* filename);
void write_xy(double* vec, double* x, int* la, char* filename);
void eig_poisson1D(double* eigval, int *la);
double eigmax_poisson1D(int *la);
double eigmin_poisson1D(int *la);
double richardson_alpha_opt(int *la);
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite);
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv);
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv);
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite);
int indexABCol(int i, int j, int *lab);
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info);
void poisson1D_CSR(int la, double *values, int *col_indices, int *row_ptr);
void poisson1D_CSC(int la, double *values, int *row_indices, int *col_ptr);
void dcscmv(int la, double *values, int *row_indices, int *col_ptr, double *x, double *y);
void dcsrmv(int la, double *values, int *col_indices, int *row_ptr, double *x, double *y);
void richardson_alpha_csc(double *values, int *row_indices, int *col_ptr, double *RHS, double *X, double *alpha_rich, int *la, double *tol, int *maxit, double *resvec, int *nbite);
int richardson_alpha_csr(const double *values, const int *col_indices, const int *row_ptr, const double *b, double *x, int n, double alpha, double tol, int maxiter);