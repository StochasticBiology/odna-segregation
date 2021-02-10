reset
#svg 480,320
unset key

set multiplot

set xrange [0.001:320]

set border 3 lw 2
set grid


set ytics nomirror
set xtics nomirror

set key bottom right

set xlabel "Time / dpc"
set ylabel "Transformed h"

# these developmental values come from the Johnston et al. eLife paper (via primary sources)
n01 = 100000.
nk1 = 300.
n02 = 300.
nk2 = 5000.
k1 = 29.
k2 = 7.
tau1 = 7./24.
tau2 = 16./24.

sstime = 40
ssn = 50
beta = 0.3*24

set xrange [0.001:320]

set border 3 lw 2
set grid

set ytics nomirror
set xtics nomirror

alpha1 = 2.*(nk1/n01)**(1./k1)
alpha2 = 2.*(nk2/n02)**(1./k2)

set key bottom right

# these functions give the variance contribution from different processes at different developmental stages
# this first function tells us about cell divisions
vh(x, alpha, n0, nc) = ((alpha+nc-1)/(alpha*n0))*((2/alpha)**(x) - 1)/((2/alpha) - 1)

# tau1 and tau2 cell divisions up to a maximum number
vh1(x) = (x <= k1*tau1 ? vh(floor(x/tau1), alpha1, n01, 1) : vh(k1, alpha1, n01, 1))
vh2(x) = (x <= k1*tau1 ? 0 : vh(floor((x-k1*tau1)/tau2), alpha2, n02, 1))

vhc1(x) = (x <= k1*tau1 ? vh(floor(x/tau1), alpha1, n01, 2) : vh(k1, alpha1, n01, 2))
vhc2(x) = (x <= k1*tau1 ? 0 : vh(floor((x-k1*tau1)/tau2), alpha2, n02, 2))

# ongoing turnover
vh3(x) = (x <= k1*tau1+k2*tau2 ? 0 : 2*beta*(x - k1*tau1-k2*tau2)/nk2)

# subsampling?
vh4(x) = (x <= sstime ? 0 : (1./ssn - 1./nk2))


transe(mu, sigma2) = log(-mu/(mu-1.)) + ((2.*mu-1.)*sigma2)/(2.*(mu-1.)**2*mu**2)
transv(mu, sigma2) = ((sigma2)/(((mu-1.)**2)*(mu**2)))

# combined model
vmodel0dev(x) = vh1(x) + vh2(x) + (x > 10 ?  2*0.000138*x : 0)
vmodel0(x) = 0.0129 + 2*0.000138*x
# single model
vmodel1(x) = 2*2.375e-4*x
# flat model
vmodel2(x) = 0.04806

unset key

set origin 0,0.
set size 0.5,1
plot (transe(0.5, vmodel0dev(x))+1.96*sqrt(transv(0.5, vmodel0dev(x)))) w filledcu y=0 lc rgbcolor "#CCCCFF", (transe(0.5, vmodel0dev(x))-1.96*sqrt(transv(0.5, vmodel0dev(x)))) w filledcu y=0 lc rgbcolor "#CCCCFF", transe(0.5, vmodel0dev(x)) lw 2 lc rgbcolor "#8888FF", "Data/hb-data.txt" u 2:4 ps 0.25 pt 7 lc rgbcolor "#000000"

set xrange [10:320]
set ylabel "V'(h)" offset 2,0,0
set y2range [0.005:0.5]
set y2label "b" offset -2,0,0
set yrange [0.005:0.5]
set logscale y2
set border 11
#set ytics ("1e-5" 1e-5, "1e-4" 1e-4, "1e-3" 1e-3, "1e-2" 1e-2, "1e-1" 1e-1, "1" 1)
#set y2tics ("1e5" 1e-5, "1e4" 1e-4, "1e3" 1e-3, "1e2" 1e-2, "1e1" 1e-1, "1" 1)
set ytics (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
set y2tics ("1000" 0.001, "500" 0.002, "200" 0.005, "100" 0.01, "50" 0.02, "20" 0.05, "10" 0.1, "5" 0.2, "2" 0.5)

set origin 0.45,0
set size 0.6,1
set format y "%.3f"
set format y2 "%.0e"

set logscale y
set xlabel "Time / dpc"
plot vmodel0(x) lw 2 lc rgbcolor "#8888FF", vmodel1(x) lw 2 lc rgbcolor "#FF8888", vmodel2(x) lw 2 lc rgbcolor "#888888", "Data/hb-age-vprime.csv" u ($1):($2) ps 0.4 pt 7 lc rgbcolor "#000000"

#> r2.log.0
#[1] 0.5872214
#> r2.log.1
#[1] 0.4317335
#> r2.log.2
#[1] -0.003370435


									     
									
