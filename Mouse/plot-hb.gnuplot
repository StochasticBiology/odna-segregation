reset
#svg 480,320
unset key

set multiplot

set size 1,0.5

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

set lmargin screen 0.15
set rmargin screen 0.8

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

# ongoing turnover using value for nu_f inferred from R code
vh5(x) = (x <= 10 ? 0 : 2*0.000138*x)

set xlabel "Time / dpc"
set ylabel "Transformed h"

transe(mu, sigma2) = log(-mu/(mu-1.)) + ((2.*mu-1.)*sigma2)/(2.*(mu-1.)**2*mu**2)
transv(mu, sigma2) = ((sigma2)/(((mu-1.)**2)*(mu**2)))

# combined model
vmodel1(x) = vh1(x)+vh2(x)+vh5(x)
# single model
vmodel0(x) = vh5(x)

unset key

set origin 0,0.5

plot (transe(0.5, vmodel1(x))+1.96*sqrt(transv(0.5, vmodel1(x)))) w filledcu y=0 lc rgbcolor "#CCCCFF", (transe(0.5, vmodel1(x))-1.96*sqrt(transv(0.5, vmodel1(x)))) w filledcu y=0 lc rgbcolor "#CCCCFF", transe(0.5, vmodel1(x)) lw 2 lc rgbcolor "#8888FF", "Data/hb-data.txt" u 2:4 ps 0.25 pt 7 lc rgbcolor "#000000"

set ylabel "V'(h)"
set y2range [1e-5:1]
set y2label "b"
set yrange [1e-5:1]
set logscale y2
set border 11
set ytics ("1e-5" 1e-5, "1e-4" 1e-4, "1e-3" 1e-3, "1e-2" 1e-2, "1e-1" 1e-1, "1" 1)
set y2tics ("1e5" 1e-5, "1e4" 1e-4, "1e3" 1e-3, "1e2" 1e-2, "1e1" 1e-1, "1" 1)
set origin 0,0
set format y "%.0e"
set format y2 "%.0e"

set logscale y
plot vmodel1(x) lw 2 lc rgbcolor "#8888FF", "Data/hb-age-vprime.csv" u ($1+21):($2) ps 0.25 pt 7 lc rgbcolor "#000000"



									     
									
