/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
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
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  int kv = *ku;
  *info = 0;

  for (int i = 0; i < *la - 1; i++) {

    if (AB[kv + i * (*lab)] == 0.0) {
      //recherche d'un nouveau pivot
      int row = i;
      double val = fabs(AB[kv + i * (*lab)]);

      if (fabs(AB[kv + (i + 1) * (*lab)]) > val) {  //compare le pivot actuel avec l'élément suivant de la ligne suivante pour trouver un meilleur pivot
        row = i + 1;
        val = fabs(AB[kv + (i + 1) * (*lab)]);
          for (int j = 0; j < *lab; j++) {  //on réalise un swap pour que l'élément suivant devienne le pivot
            double temp = AB[j + i * (*lab)];
            AB[j + i * (*lab)] = AB[j + row * (*lab)];
            AB[j + row * (*lab)] = temp;
        }
        ipiv[i] = row + 1;  //update tableau des pivots pour remplacer l'ancien pivot par le nouveau
      }
      
      if (val == 0.0) {
        *info = i + 1; //renvoie l'info de la ligne quand un pivot nul est détecté
        return *info;
      }

    }
    
    AB[kv + 1 + i * (*lab)] /= AB[kv + i * (*lab)]; //facteur multiplicateur
    AB[kv + (i + 1) * (*lab)] -= AB[kv + 1 + i * (*lab)] * AB[kv - 1 + (i + 1) * (*lab)]; //update diagonale principale
  }

  *info = 0;
  return *info;
}
