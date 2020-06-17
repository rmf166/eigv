set autoscale
unset logscale
unset label
set xtic auto
set ytic auto
set ylabel "Source and group flux" enhanced
set xlabel "Distance (cm)" enhanced
set title "Problem 1"
set size square
set xr [0:200]
set key bottom right
set style line 1 lt 1 lc rgb "blue" lw 3
set style line 2 lt 1 lc rgb "green" lw 3
plot    "p1-1-3-1.dat" using 1:2 title "lc" with lines ls 1, \
        "p1-1-3-1.dat" using 1:3 title "source" with lines ls 2
set terminal pdfcairo enhanced color dashed
set output "plot-1-1-lc.pdf"
replot
set terminal x11
