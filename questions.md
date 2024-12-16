## Exercice 3 

#### Question 1
En C comment doit on déclarer et allouer une matrice pour utiliser **BLAS** et **LAPACK** ?  
Pour utiliser **BLAS** et **LAPACK**, les matrices doivent respecter le stockage priorité colonne. Comme nous avons une matrice tridiagonale, celle-ci sera stockée dans un tableau compact, pour limiter le stockage. Ce sera alors une matrice dite matrice bande.  
L'allocation peut se faire de la manière suivante :  
```c
double* AB = (double*)malloc(M * N * sizeof(double));
```
avec M le nombre de lignes et N le nombre de colonnes.

#### Question 2
Quelle est la signification de la constante **LAPACK_COL_MAJOR** ?  
La constante **LAPACK_COL_MAJOR** fait référence au stockage en mémoire, en priorité colonne.  

#### Question 3
A quoi correspond la dimension principale généralement notée ld ?  
Lorsqu'une matrice est stockée en colonne-major, la dimension principale notée ld correspond à la distance mémoire entre deux éléments consécutifs d'une colonne. Dans le cas d'une matrice n x m, cela correspond généralement à m.  

#### Question 4
Que fait la fonction **dgbmv** ? Quelle méthode implémente-t-elle ?  
La fonction **dgbmv** dans BLAS est utilisée dans le cadre d'une multiplication matrice-vecteur d'une matrice bande. La fonction calcule alors :
$$
y = \alpha A x + \beta y
$$

#### Question 5
Que fait la fonction **dgbtrf** ? Quelle méthode implémente-t-elle ?  
La fonction **dgbtrf** de LAPACK est utilisée pour effectuer une factorisation partielle LU d'une matrice bande. La focntion calcule alors : 
$$
A = P L U
$$

#### Question 6
Que fait la fonction **dgbtrs** ? Quelle méthode implémente-t-elle ?  
La fonction **dgbtrs** de LAPACK résout le système :  
$$
A X = B 
$$
en utilisant la factorisation générée par **dgbtrf**.  

#### Question 7
Que fait la fonction **dgbsv** ? Quelle méthode implémente-t-elle ?  
La fonction dgbsv revient à faire la fonction **dgbtrf** suivie de la fonction **dgbtrs**. Elle effectue donc la factorisation LU suivie de la résolution du système en un seul appel de fonction.  

#### Question 8
Comment calculer la norme du résidu relatif avec des appels BLAS ?  
On rappelle la formule de la norme du résidu :  
$$
\text{résidu relatif} = \frac{\| A X - B \|}{\| B \|}  
$$  
On appelera alors la fonction **dgemv** pour calculer AX, puis on fait tout simplement la différence avec B.  
Ensuite, on pourra calculer la norme de AX - B, ainsi que celle de B, avec **snrm2**/**scnrm2** ou **dnrm2**/**dznrm2**.  
Le résidu se calule ensuite en faisant la division entre les deux normes obtenues.