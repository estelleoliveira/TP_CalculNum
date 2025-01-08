#include "lib_poisson1D.h"
#include <stdlib.h>
#include <stdio.h>

/*On rempli les 3 tableaux nécessaires pour le problème de poisson; les valeurs non nulles de la matrice, 
les indices de lignes correspondant et enfin le tableau de pointeur de colonnes.*/
void poisson1D_CSC(int la, double *values, int *row_indices, int *col_ptr){
    int idx = 0;    //initialisation
    col_ptr[0] = 0; //initialisation

    for (int i = 0; i < la; ++i){
        if (i > 0){
            values[idx] = -1.0;       //valeur de la diagonale inférieure
            row_indices[idx] = i - 1; //indice de la ligne correspondante
            ++idx;
        }

        values[idx] = 2.0;    //valeur de la diagonale principale
        row_indices[idx] = i; //indice de la ligne correspondante
        ++idx;

        if (i < la - 1){
            values[idx] = -1.0;       //valeur de la diagonale supérieure
            row_indices[idx] = i + 1; //indice de la ligne correspondante
            ++idx;
        }
        col_ptr[i + 1] = idx; //pointeur de colonne pour la prochaine colonne
    }
}
//produit matrice-vecteur pour CSC
void dcscmv(int la, double *values, int *row_indices, int *col_ptr, double *x, double *y){
    for (int i = 0; i < la; ++i) {
        y[i] = 0.0;
    }

    for (int i = 0; i < la; ++i){
        for (int j = col_ptr[i]; j < col_ptr[i + 1]; ++j){
            y[row_indices[j]] += values[j] * x[i];
        }
    }
}