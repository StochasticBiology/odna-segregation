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

# these values come from the R likelihood maximisation
f0(x) = 0.012977+2*x*0.0001375
f1(x) = 2*x*0.0002375
f2(x) = 0.04806

# we are now considering LE and so need the expression with selection, with parameters from R model fit
e = exp(1.)
eselect(x,rho) = 1./(1.+exp(-rho*x))
vselectfit(x, nuf, v0, rho, n) = vh1(x)+vh2(x)+(x <= 10 ? 0 : exp(-2./(1.+exp(rho*x))) * ( 4.*e*nuf + 4.*exp(1.+rho*x)*nuf + exp(2./(1.+exp(rho*x)) + rho*x)*(rho - 4.*nuf) - exp(2./(1.+exp(rho*x)))*(rho + 4.*nuf)) / (4.*(exp(rho*x)+1.)*n*rho))
vselect(x, nuf, v0, rho, n) = v0 + exp(-2./(1.+exp(rho*x))) * ( 4.*e*nuf + 4.*exp(1.+rho*x)*nuf + exp(2./(1.+exp(rho*x)) + rho*x)*(rho - 4.*nuf) - exp(2./(1.+exp(rho*x)))*(rho + 4.*nuf)) / (4.*(exp(rho*x)+1.)*n*rho)


# these values come from the R model fitting
# "all" model
emodelle3(x) = eselect(x, -0.000515423)
vmodelle3(x) = vselect(x, 0.7295, 0.00778, -0.000515423, 1517.89)
# no v0
emodelle4(x) = eselect(x, -0.00049066)
vmodelle4(x) = vselect(x, 0.7131, 0, -0.00049066, 1090.5)
# no nuf
emodelle5(x) = eselect(x, -0.00058910)
vmodelle5(x) = vselect(x, 0, 0.04152, -0.00058910, 618918.7)

transe(mu, sigma2) = log(-mu/(mu-1.)) + ((2.*mu-1.)*sigma2)/(2.*(mu-1.)**2*mu**2)
transv(mu, sigma2) = ((sigma2)/(((mu-1.)**2)*(mu**2)))

emodelhb0(x) = 0.5
vmodelhb0(x) = f0(x)
emodelhb1(x) = 0.5
vmodelhb1(x) = f1(x)
emodelhb2(x) = 0.5
vmodelhb2(x) = f2(x)

set origin 0,0.5
set key bottom left
plot (transe(emodelhb0(x), vmodelhb0(x))+1.96*sqrt(transv(emodelhb0(x), vmodelhb0(x)))) w l lw 2 lc rgbcolor "#CCCCFF" notitle, (transe(emodelhb0(x), vmodelhb0(x))-1.96*sqrt(transv(emodelhb0(x), vmodelhb0(x)))) w l lw 2 lc rgbcolor "#CCCCFF" title "(i)+(ii)+(v)",     (transe(emodelhb1(x), vmodelhb1(x))+1.96*sqrt(transv(emodelhb1(x), vmodelhb1(x)))) w l lw 2 lc rgbcolor "#FFCCCC" notitle, (transe(emodelhb1(x), vmodelhb1(x))-1.96*sqrt(transv(emodelhb1(x), vmodelhb1(x)))) w l lw 2 lc rgbcolor "#FFCCCC" title "(ii)",    (transe(emodelhb2(x), vmodelhb2(x))+1.96*sqrt(transv(emodelhb2(x), vmodelhb2(x)))) w l lw 2 lc rgbcolor "#AAAAAA" notitle, (transe(emodelhb2(x), vmodelhb2(x))-1.96*sqrt(transv(emodelhb2(x), vmodelhb2(x)))) w l lw 2 lc rgbcolor "#AAAAAA" title "(i)+(v)",    "Data/hb-data.txt" u 2:4 ps 0.25 pt 7 lc rgbcolor "#000000" title "HB data"

set origin 0.33,0.5
set key bottom left
plot (transe(emodelle3(x), vmodelle3(x))+1.96*sqrt(transv(emodelle3(x), vmodelle3(x)))) w l lw 2 lc rgbcolor "#CCCCFF" notitle, (transe(emodelle3(x), vmodelle3(x))-1.96*sqrt(transv(emodelle3(x), vmodelle3(x)))) w l lw 2 lc rgbcolor "#CCCCFF" title "(i)+(ii)+(v)",     (transe(emodelle4(x), vmodelle4(x))+1.96*sqrt(transv(emodelle4(x), vmodelle4(x)))) w l lw 2 lc rgbcolor "#FFCCCC" notitle, (transe(emodelle4(x), vmodelle4(x))-1.96*sqrt(transv(emodelle4(x), vmodelle4(x)))) w l lw 2 lc rgbcolor "#FFCCCC" title "(ii)",    (transe(emodelle5(x), vmodelle5(x))+1.96*sqrt(transv(emodelle5(x), vmodelle5(x)))) w l lw 2 lc rgbcolor "#AAAAAA" notitle, (transe(emodelle5(x), vmodelle5(x))-1.96*sqrt(transv(emodelle5(x), vmodelle5(x)))) w l lw 2 lc rgbcolor "#AAAAAA" title "(i)+(v)",    "Data/le-data.txt" u 2:4 ps 0.25 pt 7 lc rgbcolor "#000000" title "LE data"

set xrange [0:100]
set xtics (0, 10, 20, 30, 40, 50, 60, 70, "Next\ngeneration" 100) nomirror
set ytics nomirror
set yrange [0.001:0.5]
set logscale y

set key bottom right
set origin 0.66,0.5
plot  vh1(x)+vh2(x) title "(i)+(v)" lw 2 lc rgbcolor "#AAAAAA", vh1(x)+vh2(x)+vh3(x) title "(i)+(ii)+(v)" lw 2 lc rgbcolor "#FFCCCC", vhc1(x)+vhc2(x) title "(i)+(v) nc = 2" lw 2 lc rgbcolor "#99CC99", vh1(x)+vh2(x)+vh4(x) title "(i)+(iv)+(v)" lw 2 lc rgbcolor "#AA88FF",  "Data/mouse-data-elife.txt" pt 7 ps 0.5 lc rgbcolor "#000000" title "NZB-BALB/C"


unset logscale y
unset key
set style fill solid
set xrange [*:*]
set yrange [*:*]
set xtics 0.00005
set origin 0,0
set ylabel "P(nu f/n)"
set xlabel "nu f/n"
plot "nuf-hist.txt" w boxes lc rgbcolor "#CCCCFF"

unset logscale y
unset key
set style fill solid
set xtics 0.00005
set ytics 0.005
set origin 0.66,0
set xlabel "nu f/n"
set ylabel "V'0"
plot "nuf-v0-scatter.txt" pt 7 ps 0.5 lc rgbcolor "#8888FF"


set origin 0.33,0
set xtics 0.005 nomirror
set x2range [0:0.02]
set border 7 lw 2
set ytics auto
set x2label "Approx min copy number"
set x2tics ("2500" 0.005, "950" 0.01, "550" 0.015, "370" 0.02) nomirror
set ylabel "P(V'0)"
set xlabel "V'0"
plot "v0-hist.txt" w boxes lc rgbcolor "#CCCCFF"
