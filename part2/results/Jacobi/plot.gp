# Configuration générale
set terminal png enhanced size 1200,800
set grid
set key top right

# Figure 1 : Comparaison solution calculée vs solution exacte
set output 'comparison_solution_Jacobi.png'
set title 'Comparaison entre la solution calculée et la solution exacte'
set xlabel 'Points de la grille'
set ylabel 'Température'
plot 'SOL.dat' with lines lw 2 title 'Solution calculée', \
     'EX_SOL.dat' with lines lw 2 title 'Solution exacte'

# Figure 2 : Historique de convergence (résidus)
set output 'convergence_history_Jacobi.png'
set title 'Historique de convergence'
set xlabel 'Itérations'
set ylabel 'Résidu relatif'
set logscale y
plot 'RESVEC.dat' with lines lw 2 title 'Résidu'
