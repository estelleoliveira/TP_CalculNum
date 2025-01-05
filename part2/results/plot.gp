# Set the output file and the terminal type
set terminal png
set output 'jacobi_comparison.png'

# Set the title, labels, and key
set title "Convergence History Comparison"
set xlabel "Iteration"
set ylabel "Error"
set key outside

# Use logarithmic scale for both x and y axes
#set logscale x
set logscale y

# Plot data from RESVEC1.dat and RESVEC2.dat
plot 'Jacobi/RESVEC.dat' using 1 with lines title 'Jacobi method', \
     'Richardson_Alpha/RESVEC.dat' using 1 with lines title 'Richardson_Alpha method', \
     'Gauss-Seidel/RESVEC.dat' using 1 with lines title 'Gauss-Seidel method', \

# Export to png
set output