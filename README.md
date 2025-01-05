##########################################################
# README
#
# T. Dufaud
##########################################################

This directory contains the code corresponding to the solution
of Poisson 1D problem by direct method or iterative method.  
It is organized in three directories:  
src/  
include/  
bin/  

"src" contains the source codes, "include" contains the 
header files and "bin" contains the executables.  
The compilation and execution can be done using the Makefile.  

Here are the principal targets: 
testenv: bin/tp_testenv
tp2poisson1D_direct: bin/tpPoisson1D_direct
tp2poisson1D_iter: bin/tpPoisson1D_iter

The command,  
$ make target  
Compile an executable bin/target   

$ make all  
compile the executable corresponding to all targets  

$ make run_target  
Execute ./bin/target  

$ make clean  
rm *.o bin/*  

Here:  
$ make run_testenv  
$ make run_tpPoisson1D_iter  
$ make run_tpPoisson1D_direct  

Gnuplot :  
$ cd part2  
$ cd results  
$ gnuplot plot.gp   //génère un png de comparaison de convergence des méthodes utilisées  

une fois dans **results**, on peut executer les commandes suivantes  
$ cd Gauss-Seidel  
$ gnuplot plot.gp   //génère un png de convergence et de comparaison avec la solution exacte pour la méthode Gauss-Seidel  
ou  
$ cd Jacobi  
$ gnuplot plot.gp   //génère un png de convergence et de comparaison avec la solution exacte pour la méthode Jacobi  
ou  
$ cd Richardson_Alpha  
$ gnuplot plot.gp   //génère un png de convergence et de comparaison avec la solution exacte pour la méthode Richardson_Alpha  