/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

#define ALPHA 0
#define JAC 1
#define GS 2
#define CSC 3
#define CSR 4

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, lab, kv;
  int *ipiv;
  int info;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *SOL, *EX_SOL, *X;
  double *AB;
  double *MB;
  
  double temp, relres;

  double opt_alpha;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  /* Size of the problem */
  NRHS=1;
  nbpoints=12;
  la=nbpoints-2;

  /* Dirichlet Boundary conditions */
  T0=5.0;
  T1=20.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  SOL=(double *) calloc(la, sizeof(double)); 
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  //write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  
  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  /* uncomment the following to check matrix A */
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  
  /********************************************/
  /* Solution (Richardson with optimal alpha) */

  /* Computation of optimum alpha */
  opt_alpha = richardson_alpha_opt(&la);
  printf("Optimal alpha for simple Richardson iteration is : %lf",opt_alpha); 

  /* Solve */
  double tol=1e-3;
  int maxit=1000;
  double *resvec;
  int nbite=0;

  resvec=(double *) calloc(maxit, sizeof(double));

  /* Solve with Richardson alpha */
  if (IMPLEM == ALPHA) {
    richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    write_vec(resvec, &nbite, "results/Richardson_Alpha/RESVEC.dat");
    write_vec(SOL, &la, "results/Richardson_Alpha/SOL.dat");
    write_vec(EX_SOL, &la, "results/Richardson_Alpha/EX_SOL.dat");
  }

  /* Richardson General Tridiag */

  /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
  kv = 1;
  ku = 1;
  kl = 1;
  MB = (double *) malloc(sizeof(double)*(lab)*la);

  if (IMPLEM == JAC) {
    extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");
    richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    write_vec(resvec, &nbite, "results/Jacobi/RESVEC.dat");
    write_vec(SOL, &la, "results/Jacobi/SOL.dat");
    write_vec(EX_SOL, &la, "results/Jacobi/EX_SOL.dat");

  } 
  
  if (IMPLEM == GS) {
    extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");
    richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    write_vec(resvec, &nbite, "results/Gauss-Seidel/RESVEC.dat");
    write_vec(SOL, &la, "results/Gauss-Seidel/SOL.dat");
    write_vec(EX_SOL, &la, "results/Gauss-Seidel/EX_SOL.dat");
  }


  if (IMPLEM == CSC) {
    int nnz = 3 * la - 2;
    double *values = malloc(sizeof(double) * nnz);
    int *row_indices = malloc(sizeof(int) * nnz);
    int *col_ptr = malloc(sizeof(int) * (la + 1));
    double *v = calloc(la, sizeof(double));
    double *AX = calloc(la, sizeof(double));
    
    poisson1D_CSC(la, values, row_indices, col_ptr);
    
    //résidu initial: v = RHS - AX
    dcsrmv(la, values, row_indices, col_ptr, SOL, AX);
    cblas_dcopy(la, RHS, 1, v, 1);                       
    cblas_daxpy(la, -1.0, AX, 1, v, 1);
    resvec[0] = cblas_dnrm2(la, v, 1) / cblas_dnrm2(la, RHS, 1);

    double norm_res = 1.0;

    while (norm_res > tol && nbite < maxit) {
      cblas_daxpy(la, opt_alpha, v, 1, SOL, 1);
      
      dcscmv(la, values, row_indices, col_ptr, SOL, AX);
      cblas_dcopy(la, RHS, 1, v, 1);
      cblas_daxpy(la, -1.0, AX, 1, v, 1);
      
      norm_res = cblas_dnrm2(la, v, 1);
      resvec[nbite] = norm_res;
      nbite++;
    }

    write_vec(SOL, &la, "results/CSC/SOL.dat");
    write_vec(resvec, &nbite, "results/CSC/RESVEC.dat");
    write_vec(EX_SOL, &la, "results/CSC/EX_SOL.dat");

    free(values);
    free(row_indices);
    free(col_ptr);
    free(v);
    free(AX);
  }

  if (IMPLEM == CSR) {
    int nnz = 3 * la - 2;
    double *values = malloc(sizeof(double) * nnz);
    int *col_indices = malloc(sizeof(int) * nnz);
    int *row_ptr = malloc(sizeof(int) * (la + 1));
    double *v = calloc(la, sizeof(double));
    double *AX = calloc(la, sizeof(double));
    
    poisson1D_CSC(la, values, col_indices, row_ptr);
    
    //résidu initial: v = RHS - AX
    dcsrmv(la, values, col_indices, row_ptr, SOL, AX);
    cblas_dcopy(la, RHS, 1, v, 1);
    cblas_daxpy(la, -1.0, AX, 1, v, 1);
    resvec[0] = cblas_dnrm2(la, v, 1) / cblas_dnrm2(la, RHS, 1);

    double norm_res = 1.0;

    while (norm_res > tol && nbite < maxit) {
      cblas_daxpy(la, opt_alpha, v, 1, SOL, 1);
      
      dcsrmv(la, values, col_indices, row_ptr, SOL, AX);
      cblas_dcopy(la, RHS, 1, v, 1);
      cblas_daxpy(la, -1.0, AX, 1, v, 1);
      
      norm_res = cblas_dnrm2(la, v, 1);
      resvec[nbite] = norm_res;
      nbite++;
    }

    write_vec(SOL, &la, "results/CSR/SOL.dat");
    write_vec(resvec, &nbite, "results/CSR/RESVEC.dat");
    write_vec(EX_SOL, &la, "results/CSR/EX_SOL.dat");

    free(values);
    free(col_indices);
    free(row_ptr);
    free(v);
    free(AX);
  }

  /* Write solution */
  //write_vec(SOL, &la, "SOL.dat");

  /* Write convergence history */
  //write_vec(resvec, &nbite, "RESVEC.dat");

    /* Relative forward error */
  relres = relative_forward_error(SOL, EX_SOL, &la);
  printf("\nThe relative forward error is relres = %e\n",relres);


  free(resvec);
  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  free(MB);
  printf("\n\n--------- End -----------\n");
}
