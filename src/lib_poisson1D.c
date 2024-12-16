/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  //printf("la = %d\n", *la);
  //printf("lab = %d\n", *lab);
  
  for (int i = 0; i < *la * *lab; i++) {
    AB[i] = 0.0;
    //printf("kv = %d\n", *kv);
  }
  

  for (int i = 0; i < *la; i++) {
    //printf("kv : %d\n", *kv);
    AB[*kv + i * (*lab)] = -1.0;
    AB[*kv + i * (*lab) + 1] = 2.0;
    AB[*kv + i * (*lab) + 2] = -1.0; 
  }

  AB[*kv] = 0.0;
  AB[((*lab) * (*la)) - 1] = 0.0;
    
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  //printf("la = %d\n", *la);
  //printf("lab = %d\n", *lab);
  for (int i = 0; i <  *la * *lab; i++) {
    AB[i] = 0.0;
  }
  

  for (int i = 0; i < *la; i++) {
    AB[*kv + i * (*lab) + 1] = 1.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}  

void set_grid_points_1D(double* x, int* la){
}

double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  int kv = *ku + *kl;
  return *info;
}
