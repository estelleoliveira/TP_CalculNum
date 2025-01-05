#include "lib_poisson1D.h"
#include <stdlib.h>
#include <stdio.h>

//stockage CSR
/*On rempli les 3 tableaux nécessaires au stockage CSR, donc le tableau de valeurs non nulles, 
le tableau d'indice de colonnes et celui du pointeur de lignes. */
void poisson1D_CSR(int la, double *values, int *col_indices, int *row_ptr){
    int idx = 0;    //initialisation
    row_ptr[0] = 0; //intialisation

    for (int i = 0; i < la; ++i){
        if (i > 0){
            values[idx] = -1.0;       //valeur de la diagonale inférieure
            col_indices[idx] = i - 1; //indice de la colonne correspondante
            ++idx;
        }

        values[idx] = 2.0;    //valeur de la diagonale principale
        col_indices[idx] = i; //indice de la colonne correspondante
        ++idx;

        if (i < la - 1){
            values[idx] = -1.0;       //valeur de la diagonale supérieure
            col_indices[idx] = i + 1; //indice de la colonne correspondante
            ++idx;
        }
        row_ptr[i + 1] = idx; //pointeur de ligne pour la prochaine ligne
    }
}
//produit matrice-vecteur pour CSR
void dcsrmv(int la, double *values, int *col_indices, int *row_ptr, double *x, double *y){
    for (int i = 0; i < la; ++i){
        double sum = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j){
            sum += values[j] * x[col_indices[j]];
        }
        y[i] = sum;
    }
}