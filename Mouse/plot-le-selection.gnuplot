reset
# svg 480,320
unset key
set multiplot

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

set xrange [0.001:320]

set border 3 lw 2
set grid

set ytics 1 nomirror
set xtics nomirror

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

# ongoing turnover using value for nu_f inferred from R code -- note this is for HB and not used here
vh5(x) = (x <= 10 ? 0 : 2*0.000138*x)

# we are now considering LE and so need the expression with selection, with parameters from R model fit
e = exp(1.)
eselect(x,rho) = 1./(1.+exp(-rho*x))
vselectfit(x, nuf, v0, rho, n) = vh1(x)+vh2(x)+(x <= 10 ? 0 : exp(-2./(1.+exp(rho*x))) * ( 4.*e*nuf + 4.*exp(1.+rho*x)*nuf + exp(2./(1.+exp(rho*x)) + rho*x)*(rho - 4.*nuf) - exp(2./(1.+exp(rho*x)))*(rho + 4.*nuf)) / (4.*(exp(rho*x)+1.)*n*rho))
vselect(x, nuf, v0, rho, n) = v0 + exp(-2./(1.+exp(rho*x))) * ( 4.*e*nuf + 4.*exp(1.+rho*x)*nuf + exp(2./(1.+exp(rho*x)) + rho*x)*(rho - 4.*nuf) - exp(2./(1.+exp(rho*x)))*(rho + 4.*nuf)) / (4.*(exp(rho*x)+1.)*n*rho)


# these values come from the R model fitting
# "all" model
emodel3(x) = eselect(x, -0.000515423)
vmodel3(x) = vselect(x, 0.7295, 0.00778, -0.000515423, 1517.89)
vmodel3fit(x) = vselectfit(x, 0.7295, 0.00778, -0.000515423, 1517.89)
# no v0
vmodel4(x) = vselect(x, 0.7131, 0, -0.00049066, 1090.5)
# no nuf
vmodel5(x) = vselect(x, 0, 0.04152, -0.00058910, 618918.7)

set xlabel "Time / dpc"
set ylabel "Transformed h"

transe(mu, sigma2) = log(-mu/(mu-1.)) + ((2.*mu-1.)*sigma2)/(2.*(mu-1.)**2*mu**2)
transv(mu, sigma2) = ((sigma2)/(((mu-1.)**2)*(mu**2)))

vmodel1(x) = vh1(x)+vh2(x)+vh5(x)
vmodel0(x) = vh5(x)

unset key

set origin 0,0
set size 0.5,1
set yrange [-6:5]
plot (transe(emodel3(x), vmodel3fit(x))+1.96*sqrt(transv(emodel3(x), vmodel3fit(x)))) w filledcu y=0 lc rgbcolor "#CCCCFF", (transe(emodel3(x), vmodel3fit(x))-1.96*sqrt(transv(emodel3(x), vmodel3fit(x)))) w filledcu y=0 lc rgbcolor "#CCCCFF", transe(emodel3(x), vmodel3fit(x)) lw 2 lc rgbcolor "#8888FF", "Data/le-data.txt" u 2:4 ps 0.25 pt 7 lc rgbcolor "#000000"

set xrange [10:420]
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
plot vmodel3(x) lw 2 lc rgbcolor "#8888FF", vmodel4(x) lw 2 lc rgbcolor "#FF8888", vmodel5(x) lw 2 lc rgbcolor "#888888", "Data/le-age-vprime.csv" u ($1):($2) ps 0.4 pt 7 lc rgbcolor "#000000"

