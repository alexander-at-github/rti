set title "Running Time vs Thread Count"
set xlabel "Number of Threads"
set ylabel "Running Time [ns]"
set grid
set xrange [0.5:8.5]

plot \
     "long.cylinder.triangle.geometry.data" using 1:2 with lines linewidth 1 ,\
     "long.cylinder.triangle.geometry.data" using 1:2:3:4 with yerrorbars

pause -1 "Hit any key to continue"
