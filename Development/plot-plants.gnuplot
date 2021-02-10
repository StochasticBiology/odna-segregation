set logscale y
set xrange [0.3:6.1]
set yrange [10:3000]
unset key
set xtics ("Egg" 1, "Late zygote" 2, "Early embryo" 3, "Late embryo" 4, "Mature" 5)

set title "Plant development"

set xtics nomirror
set ytics nomirror
set border 3 lw 2
set grid

set xlabel "Developmental stage"
set ylabel "mtDNA copy number per cell"
plot "Data/bneck-plantj-trans.txt" u 1:2:($3-$1):($4-$2) w vectors lc rgbcolor "#4444FF", "Data/bneck-plantj.txt" u 1:3:4 w errorbars pt 7 ps 0.75 lc rgbcolor "#0000BB", "Data/bneck-plantj.txt" u ($6 == 0 ? $1+0.05 : 1/0):($3):($1 == 3 || $1 == 5 ? stringcolumn(5)." ".stringcolumn(2) : stringcolumn(5)) w labels left font "Arial,8", "Data/bneck-plantj.txt" u ($6 == 1 ? $1-0.05 : 1/0):($3):($1 == 3 || $1 == 5 ? stringcolumn(5)." ".stringcolumn(2) : stringcolumn(5)) w labels right font "Arial, 8", "Data/bneck-plantj.txt" u 1:($3/1.2):($4 == 0 ? "?" : "") w labels 
