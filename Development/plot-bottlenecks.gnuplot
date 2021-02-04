reset
set multiplot
set size 0.5,1

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

set origin 0.5,0
set xrange [0:40]
set title "Rice, Pelargonium, Pea"
unset xlabel
set xtics ("Egg" 10, "Embryo" 20, "Shoot" 30)
plot "Data/bneck-summary.csv" u ($1 == 2 ? $3 : 1/0):4 pt 7 lc rgbcolor "#0000BB"
