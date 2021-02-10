reset
#svg 640,320

set logscale y
set yrange [10:1e6]

set border 3 lw 2
set grid

set xtics nomirror
set ytics nomirror
unset key

set origin 0,0
set xrange [-2:*]
set xlabel "Days post conception"
set ylabel "mtDNA copy number per cell"
set title "Mouse germline"
plot "Data/bneck-summary.csv" u ($1 == 1 ? $3 : 1/0):4 pt 7 lc rgbcolor "#0000BB"

