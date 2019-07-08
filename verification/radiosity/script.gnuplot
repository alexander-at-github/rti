set terminal postscript eps size 3.5,2.62 enhanced color font 'Verdana, 15' linewidth 2
    #font 'Helvetica,20' linewidth 2
set output "plot.eps"

#set style data histogram
#set boxwidth 0.95

#### change default colors
#set style line 1 linecolor rgb '#8b1a0e' pointtype 1 pointsize 1 linetype 1 linewidth 2    # red
#set style line 2 linecolor rgb '#5e9c36' pointtype 6 pointsize 1 linetype 1 linewidth 2    # green
set loadpath '.'
load 'dark2.pal'
### put border more to the background
set style line 11 linecolor rgb '#808080' linetype 1
set border 3 back linestyle 11
set tics nomirror
### add grid
set style line 12 linecolor rgb '#808080' linetype 0 linewidth 1
set grid back linestyle 12

set logscale y
set format y "10^{%T}"

set title "Cylinder with Aspect Ratio 45"
set ylabel "Normalized Flux"
set xlabel "Normalized Depth"
#set grid
#set xrange [-0.4:3.6]
#set yrange [0:]
set size ratio 1

set xtics rotate by 40 right

set key font ",8"
set key box width 1.2 height 1.2

plot \
  "analytical_solution/sp.1.data" \
  using 1:2 with line linestyle 1 linewidth 1 title "Analytical 1", \
  'data/sp.1.data' \
  using 1:2 with line linestyle 2 linewidth 2 title "8d380e6, Simulated 1", \
  \
  "analytical_solution/sp.0.1.data" \
  using 1:2 with line linestyle 6 linewidth 1 title "Analytical 0.1", \
  'data/sp.0.1.v2.data' \
  using 1:2 with line linestyle 7 linewidth 2 title "Simulated 0.1"
  #'data/sp.1.v2.data' \
  #using 1:2 with line linestyle 3 linewidth 2 title "Simulated 1 v2", \
  #'data/sp.1.v3.data'  \
  #using 1:2 with line linestyle 4 linewidth 2 title "Simulated 1 v3", \
  #'data/sp.1.commit.c0a49d3.data'  \
  #using 1:2 with line linestyle 10 linewidth 1 title "c0a49d3 (pre reflection) spheres, Simulated 1", \
  #'data/sp.1.commit.c0a49d3.v2.data'  \
  #using 1:2 with line linestyle 11 linewidth 1 title "c0a49d3 (pre reflection) triangles, Simulated 1", \
  #
#  "analytical_solution/sp_0_0625.data" \
#  using 1:2 with line linestyle 4 linewidth 1 title "Analytical 0.0625", \
#  'data/sp.0.0625.data' \
#  using 1:2 with line linestyle 4 linewidth 2 title "Simulated 0.0625", \
#  "analytical_solution/sp_0_015625.data" \
#  using 1:2 with line linestyle 2 linewidth 1 title "Analytical 0.015625", \
#  'data/sp.0.015625.v4.data' \
#  using 1:2 with line linestyle 11 linewidth 2 title "Simulated v4 0.015625"
  #'data/sp.0.015625.data' \
  #using 1:2 with line linestyle 2 linewidth 2 title "Simulated 0.015625", \
  #'data/sp.0.015625.v2.data' \
  #using 1:2 with line linestyle 8 linewidth 2 title "Simulated v2 0.015625", \
  #'data/sp.0.015625.v3.data' \
  #using 1:2 with line linestyle 10 linewidth 2 title "Simulated v3 0.015625", \

  #'analytical_solution/data.data' \
  #using 1:2 with line linestyle 1 linewidth 2 title "Analytical", \
  #'' \
  #using 1:2 every 42 with points linestyle 1 pointtype 6 notitle, \

#plot  'gcc-6.3.0.d/11.06.morning.d/small.cylinder.fine.d/data.data' \
#        using (column(0)):(column(0)+1):xtic(1) linecolor rgb "white" notitle, \
#      'gcc-6.3.0.d/11.06.morning.d/small.cylinder.fine.d/data.data' \
#        using (column(0)-0.1):3 with points pointtype 7 linecolor rgb "red" title "GCC 6.3", \
#      'intel/13.06.morning.d/small.cylinder.fine.d/data.data' \
#        using (column(0)+0.1):3 with points pointtype 7 linecolor rgb "blue" title "Intel"

