#!/bin/sh
rm results.pdf
for sol in DD SC LD LC
do
  for c in 1 2 3 
  do
    rm plot.p
    echo 'set autoscale' >> plot.p
    echo 'set logscale xy' >> plot.p
    echo 'unset label' >> plot.p
    echo 'set xtic auto' >> plot.p
    echo 'set ytic auto' >> plot.p
    echo 'set grid xtic' >> plot.p
    echo 'set grid ytic' >> plot.p
    echo 'set ylabel "delta-k" enhanced' >> plot.p
    echo 'set xlabel "{/Symbol t} (mfp)" enhanced' >> plot.p
    echo 'set format y "10^{%L}"' >> plot.p
    echo 'set format x "10^{%L}"' >> plot.p
    if [ ${c} == "1" ]; then
      echo 'set title "'${sol}' (c = 0.8)"' >> plot.p
    else
      if [ ${c} == "2" ]; then
        echo 'set title "'${sol}' (c = 0.9)"' >> plot.p
      else
        echo 'set title "'${sol}' (c = 0.99)"' >> plot.p
      fi
    fi
    echo 'set size square' >> plot.p
    echo 'set pointsize 0.5' >> plot.p
    echo 'set xr [1.0e-2:10]' >> plot.p
    echo 'set yr [1.0e-10:10]' >> plot.p
    echo 'set key top left' >> plot.p
    echo 'set key width -1' >> plot.p
    echo 'set key spacing 0.75' >> plot.p
    echo 'set key font ",8"' >> plot.p
    echo 'set style line 1 lt 1 lc rgb "red" lw 1' >> plot.p
    echo 'set style line 2 lt 1 lc rgb "green" lw 1' >> plot.p
    echo 'set style line 3 lt 1 lc rgb "blue" lw 1' >> plot.p
    echo 'set style line 4 lt 1 lc rgb "magenta" lw 1' >> plot.p
    echo 'set style line 5 lt 1 lc rgb "black" lw 1' >> plot.p
    echo 'plot    "'${sol}${c}'-num.dat" using 1:2 title "N=5" with lines ls 1, \' >> plot.p
    echo '        "'${sol}${c}'-num.dat" using 1:3 title "N=10" with lines ls 2, \' >> plot.p
    echo '        "'${sol}${c}'-num.dat" using 1:4 title "N=20" with lines ls 3, \' >> plot.p
    echo '        "'${sol}${c}'-num.dat" using 1:5 title "N=40" with lines ls 4, \' >> plot.p
    echo '        "'${sol}${c}'-num.dat" using 1:5 title "N=80" with lines ls 5   ' >> plot.p
    echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
    echo 'set output "plot-'${sol}'-'${c}'.pdf"' >> plot.p
    echo 'replot' >> plot.p
    echo 'set terminal x11' >> plot.p
    gnuplot plot.p
  done
done
pdfunite plot-DD*.pdf plot-SC*.pdf plot-LD*.pdf plot-LC*.pdf results.pdf
rm plot*.pdf
