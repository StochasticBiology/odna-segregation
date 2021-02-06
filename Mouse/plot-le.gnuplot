reset
# svg 480,320
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
nuf = 0.7295
rho = -0.000515423
v0 = 0.00778
n = 1517.89
e = exp(1.)
eselect(x) = 1./(1.+exp(-rho*x))
vselect(x) = vh1(x)+vh2(x)+(x <= 10 ? 0 : exp(-2./(1.+exp(rho*x))) * ( 4.*e*nuf + 4.*exp(1.+rho*x)*nuf + exp(2./(1.+exp(rho*x)) + rho*x)*(rho - 4.*nuf) - exp(2./(1.+exp(rho*x)))*(rho + 4.*nuf)) / (4.*(exp(rho*x)+1.)*n*rho))

set xlabel "Time / dpc"
set ylabel "Transformed h"

transe(mu, sigma2) = log(-mu/(mu-1.)) + ((2.*mu-1.)*sigma2)/(2.*(mu-1.)**2*mu**2)
transv(mu, sigma2) = ((sigma2)/(((mu-1.)**2)*(mu**2)))

vmodel1(x) = vh1(x)+vh2(x)+vh5(x)
vmodel0(x) = vh5(x)

unset key

plot (transe(eselect(x), vselect(x))+1.96*sqrt(transv(eselect(x), vselect(x)))) w filledcu y=0 lc rgbcolor "#CCCCFF", (transe(eselect(x), vselect(x))-1.96*sqrt(transv(eselect(x), vselect(x)))) w filledcu y=0 lc rgbcolor "#CCCCFF", transe(eselect(x), vselect(x)) lw 2 lc rgbcolor "#8888FF", "Data/le-data.txt" u 2:4 ps 0.25 pt 7 lc rgbcolor "#000000"

