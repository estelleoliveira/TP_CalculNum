/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <string.h>
#include <time.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));

  /* LU Factorization */
  if (IMPLEM == TRF) {
    time_t starttrf, endtrf;
    starttrf = clock();
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    endtrf = clock();
    double eclapsed_time_trf = ((double) (endtrf - starttrf))/CLOCKS_PER_SEC;
    printf("Eclapsed time dgbtrf in seconds: %f\n", eclapsed_time_trf);
  }

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if (IMPLEM == TRI) {
    time_t starttridiag, endtridiag;
    starttridiag = clock();
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    endtridiag = clock();
    double eclapsed_time_tridiag = ((double) (endtridiag - starttridiag))/CLOCKS_PER_SEC;
    printf("Eclapsed time dgbtrf in seconds: %f\n", eclapsed_time_tridiag);
    
  }

  if (IMPLEM == TRI || IMPLEM == TRF){
    /* Solution (Triangular) */
    time_t startrs, endrs;
    startrs = clock();
    if (info==0){
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
      printf("\n INFO = %d\n",info);
    }
    endrs = clock();

    double eclapsed_time_rs = ((double)(endrs - startrs))/CLOCKS_PER_SEC;
    printf("Eclapsed time dgbtrs in seconds: %f\n", eclapsed_time_rs);
  }

  /* It can also be solved with dgbsv */
  if (IMPLEM == SV) {
    time_t start, end;
    start = clock();
    int lb = 2 * kl + ku + 1;
    double *ABnew = (double *) malloc(sizeof(double)* (lb * la));
    memcpy(ABnew, AB, sizeof(double) * (lb * la)); //on copie toute la matrice AB
    ipiv = (int *) calloc(la, sizeof(int));


    // Appel de dgbsv (Factorisation LU et résolution du système en une seule étape)
    dgbsv_(&la, &kl, &ku, &NRHS, ABnew, &lb, ipiv, RHS, &la, &info);
    end = clock();

    if (info!=0){printf("\n INFO DGBSV = %d\n",info);}
    else{
      printf("\n INFO = %d\n",info);
    }
    double eclapsed_time = ((double)(end - start))/CLOCKS_PER_SEC;
    printf("Eclapsed time dgbsv in seconds: %f\n", eclapsed_time);
    free(ipiv);
    free(ABnew);
  }

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
