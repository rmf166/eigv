set autoscale
set logscale xy
unset label
set xtic auto
set ytic auto
set grid xtic
set grid ytic
set ylabel "delta-k" enhanced
set xlabel "{/Symbol t} (mfp)" enhanced
set format y "10^{%L}"
set format x "10^{%L}"
set title "LC (c = 0.99)"
set size square
set pointsize 0.5
set xr [1.0e-2:10]
set yr [1.0e-10:10]
set key top left
set key width -1
set key spacing 0.75
set key font ",8"
set style line 1 lt 1 lc rgb "red" lw 1
set style line 2 lt 1 lc rgb "green" lw 1
set style line 3 lt 1 lc rgb "blue" lw 1
set style line 4 lt 1 lc rgb "magenta" lw 1
set style line 5 lt 1 lc rgb "black" lw 1
plot    "LC3-num.dat" using 1:2 title "N=5" with lines ls 1, \
        "LC3-num.dat" using 1:3 title "N=10" with lines ls 2, \
        "LC3-num.dat" using 1:4 title "N=20" with lines ls 3, \
        "LC3-num.dat" using 1:5 title "N=40" with lines ls 4, \
        "LC3-num.dat" using 1:5 title "N=80" with lines ls 5   
set terminal pdfcairo enhanced color dashed
set output "plot-LC-3.pdf"
replot
set terminal x11
