/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>

void eig_poisson1D(double* eigval, int *la){
  double h = 1.0 / (*la + 1);  //pas entre les points de la grille
  eigval[*la];
  for (int k = 0; k < *la; k++){
    eigval[k + 1] = 4 * (sin(((double) k + 1.0) * M_PI * h / 2))*(sin(((double) k + 1.0) * M_PI * h / 2));
    printf("eigval[%d] = %f\n", k, eigval[k]);
  }
  //*eigval = 2.0/(eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

double eigmax_poisson1D(int *la){
  double h = 1.0 / (*la + 1);
  return 4 * (sin(*la * M_PI * h / 2))*(sin(*la * M_PI * h / 2)); //formule pour déterminer la valeur propre maximale dans le cas de matrice symétrique
}

double eigmin_poisson1D(int *la){
  double h = 1.0 / (*la + 1);
  return 4 * (sin(M_PI * h / 2))*(sin(M_PI * h / 2)); //formule pour déterminer la valeur propre minimale dans le cas de matrice symétrique
}

double richardson_alpha_opt(int *la){
  return 2.0/(eigmin_poisson1D(la) + eigmax_poisson1D(la)); //formule pour déterminer alpha optimal grâce aux valeurs propres min et max déterminées plus tôt
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
 double* v = (double*)calloc(*la, sizeof(double));  //on génère un vecteur v 
    if (!v) {
      fprintf(stderr, "Memory allocation failed in richardson_alpha\n");  //erreur si l'allocation n'a pas fonctionné
      return;
    }

  //calcul du résidu
  const double norm = 1.0 / cblas_dnrm2(*la, RHS, 1); // 1/norm(RHS)
  cblas_dcopy(*la, RHS, 1, v, 1); //on copie le vecteur RHS dans v
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, v, 1);  //on fait dgbmv, vi = vi - (AB*X), ce qui permet de stocker le résidu dans v; vi dans le terme de droite est RHS
  
  resvec[*nbite] = cblas_dnrm2(*la, v, 1) * norm; //on fait donc norm(v)/norm(RHS), on obtient une norme relative qu'on stocke dans resvec; donc norm(résidu)/norm(RHS) le résidu relatif
  
  while (resvec[*nbite] > *tol && *nbite < *maxit) {
    cblas_daxpy(*la, *alpha_rich, v, 1, X, 1);  //alpha_rich est utilisé pour calculer le daxpy; X = X + 0.5 * v
    
    //maj du résidu
    cblas_dcopy(*la, RHS, 1, v, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, v, 1);
    
    (*nbite)++;
    resvec[*nbite] = cblas_dnrm2(*la, v, 1) * norm;
  }

  free(v);

}

/*Permet d'extraire des diagonales de AB pour les stocker dans la matrice MB qui est la matrice préconditionné pour la méthode de Jacobi*/
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  //M = D
  for (int i = 0; i < *la; i++) {
    MB[i * (*lab) + 1] = AB[i * (*lab) + 1];  //on stocke seulement la diagonale principale
  }
}

/*Permet d'extraire des diagonales de AB pour les stocker dans la matrice MB qui est la matrice préconditionné pour la méthode de Gauss Seidel*/
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  //M = D - E
  for (int i = 0; i < *la; i++) {
    MB[i * (*lab) + 1] = AB[i * (*lab) + 1];  //on stocke la diagonale principale dans MB
    MB[i * (*lab) + 2] = AB[i * (*lab) + 2];  //on stocke la sur-diagonale dans MB
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double* v = (double*)malloc((*la) * sizeof(double));
  double* AX = (double*)malloc((*la) * sizeof(double));
  if (!v || !AX) {
    free(v);
    free(AX);
    fprintf(stderr, "Memory allocation failed in richardson_MB\n"); //erreur si v ou AX n'a pas réussi à être alloué
    return;
  }

  //initialise X à 0
  memset(X, 0, (*la) * sizeof(double));
  
  double norm_res = 1.0;
  
  while (norm_res > *tol && *nbite < *maxit) {
    //calcul résidu R = RHS - AB*X
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, AX, 1);
    cblas_dcopy(*la, RHS, 1, v, 1);
    cblas_daxpy(*la, -1.0, AX, 1, v, 1);
    
    //maj Xi = Xi * vi / MBi
    for (int i = 0; i < *la; i++) {
      X[i] += v[i] / MB[i * (*lab) + *kl];  //on utilise la matrice préconditionnée de Jacobi ou Gauss-Seidel
    }

    //calcul du nouveau résidu
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, AX, 1);
    cblas_dcopy(*la, RHS, 1, v, 1);
    cblas_daxpy(*la, -1.0, AX, 1, v, 1);
    
    norm_res = cblas_dnrm2(*la, v, 1);
    resvec[*nbite] = norm_res;
    (*nbite)++;
  }

  free(v);
  free(AX);
}

