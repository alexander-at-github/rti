set title "Running Time vs Data Structures on a mesh with 7145 vertices and 14286 triangles"
set xlabel "Trinagles / Discs / Spheres"
set ylabel "Running Time [s]"
#set grid
set xrange [-0.5:2.5]
set yrange [0:]

plot 'long.cylinder.tri.disc.sphere.data' using 0:3:4:5:xtic(1) with yerrorbars linewidth 2
#     "long.cylinder.triangle.geometry.data" using 1:2 with lines linewidth 1 ,\
#     "long.cylinder.triangle.geometry.data" using 1:2:3:4 with yerrorbars

pause -1 "Hit any key to continue"
