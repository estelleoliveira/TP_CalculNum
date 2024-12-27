/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
  double h = 1.0 / (*la + 1);  //pas entre les points de la grille
  eigval[*la];
  for (int k = 0; k < *la; k++){
    eigval[k + 1] = 4 * (sin(((double) k + 1.0) * M_PI * h / 2))*(sin(((double) k + 1.0) * M_PI * h / 2));
    printf("eigval %f\n", eigval[k]);
  }
}

double eigmax_poisson1D(int *la){
  double h = 1.0 / (*la + 1);
  return 4 * (sin(*la * M_PI * h / 2))*(sin(*la * M_PI * h / 2));
}

double eigmin_poisson1D(int *la){
  double h = 1.0 / (*la + 1);
  return 4 * (sin(M_PI * h / 2))*(sin(M_PI * h / 2));
}

double richardson_alpha_opt(int *la){
  return 2.0/(eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

