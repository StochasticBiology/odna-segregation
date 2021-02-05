reset
#svg 1024,800
unset key

# these developmental values come from the Johnston et al. eLife paper (via primary sources)
n01 = 100000.
nk1 = 500.
n02 = 500.
nk2 = 50000.
k1 = 29.
k2 = 7.
tau1 = 7./24.
tau2 = 16./24.

sstime = 40
ssn = 50
beta = 0.3*24

set xrange [0:105]

set border 3 lw 2
set grid
set xtics nomirror
set ytics nomirror

#set logscale y

alpha1 = 2.*(nk1/n01)**(1./k1)
alpha2 = 2.*(nk2/n02)**(1./k2)

set key bottom right

# these functions give the variance contribution from different processes at different developmental stages
# this first function tells us about cell divisions
vh(x, alpha, n0, nc) = ((alpha+nc-1)/(alpha*n0))*((2/alpha)**(x+1) - 1)/((2/alpha) - 1)

# tau1 and tau2 cell divisions up to a maximum number
vh1(x) = (x <= k1*tau1 ? vh(floor(x/tau1), alpha1, n01, 1) : vh(k1, alpha1, n01, 1))
vh2(x) = (x <= k1*tau1 ? 0 : vh(floor((x-k1*tau1)/tau2), alpha2, n02, 1))

vhc1(x) = (x <= k1*tau1 ? vh(floor(x/tau1), alpha1, n01, 2) : vh(k1, alpha1, n01, 2))
vhc2(x) = (x <= k1*tau1 ? 0 : vh(floor((x-k1*tau1)/tau2), alpha2, n02, 2))

# ongoing turnover
vh3(x) = (x <= k1*tau1+k2*tau2 ? 0 : 2*beta*(x - k1*tau1-k2*tau2)/nk2)

# subsampling?
vh4(x) = (x <= sstime ? 0 : (1./ssn - 1./nk2))
									

set multiplot
set size 0.33,0.5

set key bottom left
set xrange [20:320]
set xlabel "Time / dpc"
set ylabel "Transformed h"

unset key

set origin 0,0.5
#plot "hb-data.txt" u 2:4 title "HB" ps 0.3 pt 7 lc rgbcolor "#000000", "model-neutral-0.txt" u 1:2 w l lw 3 lc rgbcolor "#FF0000" notitle, "" u 1:3 w l lw 3 lc rgbcolor "#FF8888" title "(i)+(ii)+(v)", "" u 1:4 w l lw 3 lc rgbcolor "#FF8888" notitle,    "model-neutral-1.txt" u 1:2 w l lw 3 lc rgbcolor "#0000FF" notitle , "" u 1:3 w l lw 3 lc rgbcolor "#8888FF" title "(ii)", "" u 1:4 w l lw 3 lc rgbcolor "#8888FF" notitle,    "model-neutral-2.txt" u 1:2 w l lw 3 lc rgbcolor "#000000" notitle, "" u 1:3 w l lw 3 lc rgbcolor "#AAAAAA" title "(i)+(v)", "" u 1:4 w l lw 3 lc rgbcolor "#AAAAAA" notitle
plot "Data/hb-data.txt" u 2:4 title "HB" ps 0.3 pt 7 lc rgbcolor "#000000", "models-neutral.txt" u 1:2 w l lw 3 lc rgbcolor "#FF0000" notitle, "" u 1:($2+1.96*$3) w l lw 3 lc rgbcolor "#FF8888" title "(i)+(ii)+(v)", "" u 1:($2-1.96*$3) w l lw 3 lc rgbcolor "#FF8888" notitle,    "" u 1:4 w l lw 3 lc rgbcolor "#0000FF" notitle , "" u 1:($4+1.96*$5) w l lw 3 lc rgbcolor "#8888FF" title "(ii)", "" u 1:($4-1.96*$5) w l lw 3 lc rgbcolor "#8888FF" notitle,    "" u 1:6 w l lw 3 lc rgbcolor "#000000" notitle, "" u 1:($6+1.96*$7) w l lw 3 lc rgbcolor "#AAAAAA" title "(i)+(v)", "" u 1:($6-1.96*$7) w l lw 3 lc rgbcolor "#AAAAAA" notitle

# these values come from the R likelihood maximisation
f0(x) = 0.012977+2*x*0.0001375
f1(x) = 2*x*0.0002375
f2(x) = 0.04806

set logscale y
set yrange [1e-3:1]
set origin 0.33,0.5
set key bottom right
set ylabel "V'(h)"
plot "Data/hb-age-vprime.csv" u ($1+21):($2) ps 0.5 pt 7 lc rgbcolor "#000000" title "HB", "models-neutral.txt" u 1:(f0($1)) w l lw 3 lc rgbcolor "#FF8888" title "(i)+(ii)+(v)", "" u 1:(f1($1)) w l lw 3 lc rgbcolor "#8888FF" title "(ii)", "" u 1:(f2($1)) w l lw 3 lc rgbcolor "#AAAAAA" title "(i)+(v)"

set xrange [0:100]
set xtics (0, 10, 20, 30, 40, 50, 60, 70, "Next\ngeneration" 100) nomirror
set ytics nomirror

set origin 0.66,0.5
plot  "Data/mouse-data-elife.txt" pt 7 ps 0.5 lc rgbcolor "#000000" title "NZB-BALB/C", vh1(x)+vh2(x) title "(i)+(v)" lw 3 lc rgbcolor "#AAAAAA", vh1(x)+vh2(x)+vh3(x) title "(i)+(ii)+(v)" lw 3 lc rgbcolor "#FF8888", vhc1(x)+vhc2(x) title "(i)+(v) nc = 2" lw 3 lc rgbcolor "#FF8800", vh1(x)+vh2(x)+vh4(x) title "(i)+(iv)+(v)" lw 3 lc rgbcolor "#8800FF"


unset logscale y
unset key
set style fill solid
set xrange [*:*]
set yrange [*:*]
set xtics 0.00005
set origin 0,0
set ylabel "P(nu f/n)"
set xlabel "nu f/n"
plot "nuf-hist.txt" w boxes

unset logscale y
unset key
set style fill solid
set xtics 0.00005
set ytics 0.005
set origin 0.66,0
set xlabel "nu f/n"
set ylabel "V'0"
plot "nuf-v0-scatter.txt" 


set origin 0.33,0
set xtics 0.005 nomirror
set x2range [0:0.02]
set ytics auto
set x2label "Approx min copy number"
set x2tics ("2500" 0.005, "950" 0.01, "550" 0.015, "370" 0.02) nomirror
set ylabel "P(V'0)"
set xlabel "V'0"
plot "v0-hist.txt" w boxes
