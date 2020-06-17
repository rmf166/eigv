#!/bin/sh
rm results.pdf
for src in 1 1
do
  for sol in dd sc ld lc
  do
    for prb in 1
    do
      rm plot.p
      echo 'set autoscale' >> plot.p
      echo 'unset logscale' >> plot.p
      echo 'unset label' >> plot.p
      echo 'set xtic auto' >> plot.p
      echo 'set ytic auto' >> plot.p
      echo 'set ylabel "Source and group flux" enhanced' >> plot.p
      echo 'set xlabel "Distance (cm)" enhanced' >> plot.p
      echo 'set title "Problem '${prb}'"' >> plot.p
      echo 'set size square' >> plot.p
      echo 'set xr [0:200]' >> plot.p
      echo 'set key bottom right' >> plot.p
      echo 'set style line 1 lt 1 lc rgb "blue" lw 3' >> plot.p
      echo 'set style line 2 lt 1 lc rgb "green" lw 3' >> plot.p
      if [ "${sol}" = "dd" ]
      then
      echo 'plot    "p'${prb}'-'${src}'-0-1.dat" using 1:2 title "'${sol}'" with lines ls 1, \' >> plot.p
      echo '        "p'${prb}'-'${src}'-0-1.dat" using 1:3 title "source" with lines ls 2' >> plot.p
      fi
      if [ "${sol}" = "sc" ]
      then
      echo 'plot    "p'${prb}'-'${src}'-1-1.dat" using 1:2 title "'${sol}'" with lines ls 1, \' >> plot.p
      echo '        "p'${prb}'-'${src}'-1-1.dat" using 1:3 title "source" with lines ls 2' >> plot.p
      fi
      if [ "${sol}" = "ld" ]
      then
      echo 'plot    "p'${prb}'-'${src}'-2-1.dat" using 1:2 title "'${sol}'" with lines ls 1, \' >> plot.p
      echo '        "p'${prb}'-'${src}'-2-1.dat" using 1:3 title "source" with lines ls 2' >> plot.p
      fi
      if [ "${sol}" = "lc" ]
      then
      echo 'plot    "p'${prb}'-'${src}'-3-1.dat" using 1:2 title "'${sol}'" with lines ls 1, \' >> plot.p
      echo '        "p'${prb}'-'${src}'-3-1.dat" using 1:3 title "source" with lines ls 2' >> plot.p
      fi
      echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
      echo 'set output "plot-'${src}'-'${prb}'-'${sol}'.pdf"' >> plot.p
      echo 'replot' >> plot.p
      echo 'set terminal x11' >> plot.p
      gnuplot plot.p
    done  
  done
done
pdfunite plot*.pdf results.pdf
rm plot*.pdf
