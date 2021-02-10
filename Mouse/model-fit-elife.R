# these developmental values come from the Johnston et al. eLife paper (via primary sources)
n01 = 100000.
nk1 = 500.
n02 = 500.
nk2 = 50000.
k1 = 29.
k2 = 7.
tau1 = 7./24.
tau2 = 16./24.

alpha1 = 2.*(nk1/n01)**(1./k1)
alpha2 = 2.*(nk2/n02)**(1./k2)

# different functions generating V'(h) predictions from theory
vh <- function(x, alpha, n0, nc) {
  return ( ((alpha+nc-1)/(alpha*n0))*((2/alpha)**(x+1) - 1)/((2/alpha) - 1) )
}

# tau1 and tau2-length cell cycles up to a maximum number k1, k2
vh1 <- function(x,sstime,ssn,beta,nc) {
  if(x <= k1*tau1) { return(vh(floor(x/tau1), alpha1, n01, 1)) }
  else { return(vh(k1, alpha1, n01, 1)) }
}

vh2 <- function(x,sstime,ssn,beta,nc) {
  if(x <= k1*tau1) { return(0) }
  else { return(vh(floor((x-k1*tau1)/tau2), alpha2, n02, 1)) }
}

# inheritance of clusters of size nc
vhc1 <- function(x,sstime,ssn,beta,nc) {
  if(x <= k1*tau1) { return(vh(floor(x/tau1), alpha1, n01, nc)) }
  else { return(vh(k1, alpha1, n01, nc)) }
}

vhc2 <- function(x,sstime,ssn,beta,nc) {
  if(x <= k1*tau1) { return(0) }
  else { return(vh(floor((x-k1*tau1)/tau2), alpha2, n02, nc)) }
}

# ongoing turnover with rate beta
vh3 <- function(x,sstime,ssn,beta,nc) {
  if(x <= k1*tau1+k2*tau2) { return(0) }
  else { return(2*beta*(x - k1*tau1-k2*tau2)/nk2) }
}

# subsampling to ssn at sstime
vh4 <- function(x,sstime,ssn,beta,nc) {
  if(x <= sstime) { return(0) }
  else { return((1./ssn - 1./nk2)) }
}

# overall model "likelihood" function -- really just a sum of squared differences
model.lik <- function(t,Y,theta,model,output=0) {
  # interpret parameters (exps prevent negative values)
  sstime.p = exp(theta[1])
  ssn.p = exp(theta[2])
  beta.p = theta[3]
  nc.p = exp(theta[4])

  # compute contributions from each mechanism
  p1 = sapply(t, vh1, sstime=sstime.p, ssn=ssn.p, beta=beta.p, nc=1)
  p2 = sapply(t, vh2, sstime=sstime.p, ssn=ssn.p, beta=beta.p, nc=1)
  p1c = sapply(t, vhc1, sstime=sstime.p, ssn=ssn.p, beta=beta.p, nc=nc.p)
  p2c = sapply(t, vhc2, sstime=sstime.p, ssn=ssn.p, beta=beta.p, nc=nc.p)
  p3 = sapply(t, vh3, sstime=sstime.p, ssn=ssn.p, beta=beta.p, nc=nc.p)
  p4 = sapply(t, vh4, sstime=sstime.p, ssn=ssn.p, beta=beta.p, nc=nc.p)

  # output either sum of squares or lists of predicted values
  # structure determined by model choice
  if(output == 0) {
    if(model == 0) { return(sum((log(p1+p2)-log(Y))**2)) }
    if(model == 1) { return(sum((log(p1+p2+p3)-log(Y))**2)) }
    if(model == 2) { return(sum((log(p1c+p2c)-log(Y))**2)) }
    if(model == 3) { return(sum((log(p1+p2+p4)-log(Y))**2)) }
  }
  else {
    if(model == 0) { return(p1+p2) }
    if(model == 1) { return(p1+p2+p3) }
    if(model == 2) { return(p1c+p2c) }
    if(model == 3) { return(p1+p2+p4) }
  }
}

# get data and propose starting params
df = read.table("Data/mouse-data-elife.txt")
ic.theta = c(log(60), log(50), 0.3*24, log(1))
model.lik(df[,1], df[,2], ic.theta, 3)

# fit model variants
opt.model.0 = optim(ic.theta, model.lik, t=df[,1], Y=df[,2], model=0, output=0)
opt.model.1 = optim(ic.theta, model.lik, t=df[,1], Y=df[,2], model=1, output=0)
opt.model.2 = optim(ic.theta, model.lik, t=df[,1], Y=df[,2], model=2, output=0)
opt.model.3 = optim(ic.theta, model.lik, t=df[,1], Y=df[,2], model=3, output=0)

# get R2s
r2.log.0 = 1 - opt.model.0$value / sum((log(df[,2])-mean(log(df[,2])))**2)
r2.log.1 = 1 - opt.model.1$value / sum((log(df[,2])-mean(log(df[,2])))**2)
r2.log.2 = 1 - opt.model.2$value / sum((log(df[,2])-mean(log(df[,2])))**2)
r2.log.3 = 1 - opt.model.3$value / sum((log(df[,2])-mean(log(df[,2])))**2)

# get best parameters
theta.0 = opt.model.0$par
theta.1 = opt.model.1$par
theta.2 = opt.model.2$par
theta.3 = opt.model.3$par

# sanity check plot
t.seq = seq(0, 100)
line.0 = model.lik(t.seq, 0, theta.0, model=0,output=1)
line.1 = model.lik(t.seq, 0, theta.1, model=1,output=1)
line.2 = model.lik(t.seq, 0, theta.2, model=2,output=1)
line.3 = model.lik(t.seq, 0, theta.3, model=3,output=1)

plot(df[,1], log(df[,2]))
points(t.seq, log(line.0), type="l")
points(t.seq, log(line.1), type="l")
points(t.seq, log(line.2), type="l")
points(t.seq, log(line.3), type="l")
