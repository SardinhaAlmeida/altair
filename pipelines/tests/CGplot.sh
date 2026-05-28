#!/bin/sh
echo 'reset
set terminal pdfcairo enhanced font "Verdana,12"
set output "CGplot.pdf"
set style line 11 lc rgb "#808080" lt 1
set border 3 back ls 11
set tics nomirror
set size ratio 1
set style line 12 lc rgb "#808080" lt 0 lw 2
set grid back ls 12
set style line 1 lc rgb "#8b1a0e" pt 1 ps 1 lt 1 lw 3
set style line 2 lc rgb "#5e9c36" pt 6 ps 1 lt 1 lw 3
set key bottom right
set xlabel "k-mer"
set ylabel "Distribution (%)"
set xtics 1
set xrange [10:13]
set yrange [0:100]
plot "CG-data.eg" u 1:2 t "AT" w lp ls 1,  "CG-data.eg" u 1:3 t "CG" w lp ls 2 ' | gnuplot -persist
