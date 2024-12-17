/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){

  for (int i = 0; i < *la; i++) {
    if (*kv > 0) {
      AB[i * (*lab)] = 0.0; //si kv = 0, on ne rajoute pas de ligne de 0.0 à notre matrice bande
    }
    AB[*kv + i * (*lab)] = -1.0;    //correspond à A[0] pour kv = 0 et A[1] pour kv = 1 (i=0)
    AB[*kv + i * (*lab) + 1] = 2.0;
    AB[*kv + i * (*lab) + 2] = -1.0; 
    
  }

  AB[*kv] = 0.0;  //remplace la première valeur par 0.0 (ici on met kv pour remplacer la valeur 0 ou 1)
  AB[((*lab) * (*la)) - 1] = 0.0; //remplace la dernière valeur par 0.0
  
  /*for (int j = 0; j < *la * *lab; j++) {
    printf("AB : %f\n", AB[j]);
    printf("kv : %d\n", *kv);
  }*/
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  
//on utilise le même fonctionnement que pour la fonction précedente en modifiant les valeurs remplies
  for (int i = 0; i < *la; i++) {
    if (*kv > 0) {
      AB[i * (*lab)] = 0.0;
    }
    AB[*kv + i * (*lab)] = 0.0;
    AB[*kv + i * (*lab) + 1] = 1.0;
    AB[*kv + i * (*lab) + 2] = 0.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  RHS[0]= *BC0; //température condition de bord à x = 0
  RHS[*la - 1]= *BC1; //température condition de bord à x = 1
  for (int i = 1; i < *la - 1; i++){
    RHS[i]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  for (int i = 0; i < *la; i++) {
    EX_SOL[i] = *BC0 + (*BC1 - *BC0) * X[i];  //T = T0 + (T1 - T0) * X
  }
}  

void set_grid_points_1D(double* x, int* la){
  double h = 1.0 / (*la + 1);  //pas entre les points de la grille
  for (int i = 0; i < *la; i++) {
    x[i] = (i + 1) * h;
  }
}

double relative_forward_error(double* x, double* y, int* la){

  double normdiff = 0.0;
  double norm_y = 0.0;
  //on effectue l'erreur relative de la forme [norm2(x - y)/norm2(y)] où x = AX et y = B
  for (int i = 0; i < *la; i++) {
    double diff = x[i] - y[i];
    normdiff += diff * diff;
    //norm_y += y[i] * y[i];
  }

  norm_y = cblas_dnrm2(*la, y, 1);
  if (norm_y == 0.0) {
    return -1.0;  //permet de vérifier que on ne va pas diviser par zéro
  }
  double erreur_relative = sqrt(normdiff) / norm_y;
  
  
  return erreur_relative;
}

int indexABCol(int i, int j, int *lab){
  return i + j * *lab;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  int kv = *ku + *kl;
  return *info;
}
