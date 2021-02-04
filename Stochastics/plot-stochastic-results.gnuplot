reset
set multiplot

# svg 600,360
set style fill solid
set rmargin screen 0.8
set bmargin screen 0.2
set tmargin screen 0.95
unset key
set size 0.8,1
set boxwidth 0.1
set xrange [-0.35:8]
set yrange [-0.0005:*]
set border 3 lw 2
set grid
set xtics nomirror
set ytics nomirror
set ylabel "V'(h)"
set xtics ("No random\ninfluence" 0, "Division\nnc = 1" 1, "Division\nnc = 1\n+ reamp" 2, "Division\nnc = 4\n+ reamp" 3, "Subsample\n25%%\n+reamp" 4, "High fusion\n(recomb)" 5, "High fission\n(turnover)" 6, "Balanced\n(both)" 7) rotate 

myrad = 0.1
nexpt = 9
set label 1 at 0, 0.0005 "N = 500" rotate font "Arial,9" tc rgbcolor "#000000"
set label 2 at 0.2, 0.0005 "1000" rotate font "Arial,9" tc rgbcolor "#000000"
set label 3 at 0.4, 0.0005 "1500" rotate font "Arial,9" tc rgbcolor "#000000"
set label 4 at 0.6, 0.0005 "2000" rotate font "Arial,9" tc rgbcolor "#000000"

set yrange [-0.0005:0.014]
set ytics (0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014)
set y2range [-0.0005:0.014]
set y2tics ("∞" 0, "500" 0.002, "250" 0.004, "167" 0.006, "125" 0.008, "100" 0.01, "83" 0.012, "71" 0.014)
set border 11
set y2label "b"

# the columns in the output file are
# 1 expt ref ; 2 popn ; 3 t ; 4 meanh ; 5 varh ; 6 var'h ; 7 meann ; 8 varn ; 9 meanf ; 10 varf ; 11 predicted mean h ; 12 predicted var h ; 13 predicted var' h
# to look at dV'(h) / dt we plot $1 (plus offset from $2) against $13 with boxes, and against $6 with circles
plot "stochastic-simulation-output.txt" u (($1 <= 4 && $3 == 1 ? $1-0.1+$2/3000 : 1/0)):13 w boxes lc rgbcolor "#888888", "" u (($1 <= 4 && $3 == 1 ? $1-0.1+$2/3000 : 1/0)):6:(myrad) w circles lc rgb "#FFFFFF" fill solid border lt 0 notitle, "" u (($1 > 4 && $3 == 1 ? $1-0.1+$2/3000 : 1/0)):13 w boxes lc rgbcolor "#000000", "" u (($1 > 4 && $3 == 1? $1-0.1+$2/3000 : 1/0)):6:(myrad) w circles lc rgb "#FFFFFF" fill solid border lt 0 notitle

unset xtics
unset ytics
unset y2tics
unset ylabel
set y2label "P(h)"

## next, a couple of insets demonstrating distribution widths for different V'(h)
set yrange [0:*]
set xrange [0:1]
set xtics (0, 1) nomirror
set lmargin screen 0.87
set tmargin screen 0.95
set rmargin screen 0.95
set bmargin screen 0.8

unset label 1
unset label 2
unset label 3
unset label 4
sigma = 0.015
set border 9
plot 1./sqrt(2.*3.1415*sigma/4)*exp(-(x-0.5)**2/(2*sigma/4)) w filledcu lc rgbcolor "#AAAAAA"

set origin 0.75,0.2
set lmargin screen 0.87
set bmargin screen 0.225
set rmargin screen 0.95
set tmargin screen 0.375

sigma = 0.0001
plot 1./sqrt(2.*3.1415*sigma/4)*exp(-(x-0.5)**2/(2*sigma/4)) w filledcu lc rgbcolor "#AAAAAA"

## finally, an inset showing the comparison between gene conversion and turnover contributions as f changes

set lmargin screen 0.6
set bmargin screen 0.7
set rmargin screen 0.73
set tmargin screen 0.9

set border 3

nu = 1
kappa = 0.002
set ytics 0.002 nomirror
set xlabel "f"
unset y2label
set ylabel "V'(h)" offset 1,0
set label 1 at  0.125, 0.0045 "R" font "Arial,9" tc rgbcolor "#888888"
set label 2 at  0.47, 0.0045 "T (N = 500)" font "Arial,9" tc rgbcolor "#000000"
set label 3 at  0.7, 0.0015 "N = 2000" font "Arial,9" tc rgbcolor "#000000"
plot 2*nu*x/500 lc rgbcolor "#000000" lw 2, 2*nu*x/2000 lc rgbcolor "#000000" lw 2, 2*kappa*(1-x)**2 lc rgbcolor "#AAAAAA" lw 2
