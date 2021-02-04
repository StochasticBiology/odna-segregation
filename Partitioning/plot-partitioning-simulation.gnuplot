reset
set multiplot
set size 0.5

# svg 640,640

set xlabel "N"

set logscale y
set format y "%.0e" 

set origin 0,0
set ylabel "V(N)"
plot "partitioning-stats-out.txt" u ($2 == 0 ? $1 : 1/0):($4*$4) pt 1 lc rgbcolor "#FF8888" title "Random", "" u ($2 == 1 ? $1 : 1/0):($4*$4) pt 2 lc rgbcolor "#8888FF" title "Spaced", "" u 1:(0.25/$1) w l lc rgbcolor "#888888" title "Binomial" #, "" u 1:(0.1/($1**1.25)) w l lc rgbcolor "#888888" title "Power-law estimate"

set origin 0.5,0.
set ylabel "V'(h)"
plot "partitioning-stats-out.txt" u ($2 == 0 ? $1 : 1/0):($6*$6/0.25) pt 1 lc rgbcolor "#FF8888" title "Random", "" u ($2 == 1 ? $1 : 1/0):($6*$6/0.25) pt 2 lc rgbcolor "#8888FF" title "Spaced", "" u 1:(1./$1) w l lc rgbcolor "#888888"  title "Hypergeometric"

unset logscale y
set origin 0,0.5
set border 0
unset xtics
unset ytics
unset ztics
unset xlabel
unset ylabel
splot "ex-389-0.txt" pt 7 ps 1 lc rgbcolor "#AAFF8888" title "Random"

set origin 0.5,0.5
splot "ex-389-1.txt" pt 7 ps 1 lc rgbcolor "#AA8888FF" title "Spaced"
