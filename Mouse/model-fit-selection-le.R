model.lik <- function(theta, t, Y, model = 0, output=1) {

  # grab the model variables from the parameter vector
  # we use logged values to ensure they all remain non-negative
  nuf = exp(theta[1])
  v0 = exp(theta[2])
  rho = theta[3]
  n = exp(theta[4])

  if(model < 3) {
    mu = 1/2+(t-t)
    if(model == 1) { v0 = 0 }
    if(model == 2) { nuf = 0 }
    sigma = sqrt(v0 + nuf*2*t)
  }
  else {
    if(model == 4) { v0 = 0 }
    if(model == 5) { nuf = 0 }
    mu = 1/(1 + exp(-rho*t))
    sigma = sqrt( v0 + (exp(-2/(1+exp(rho*t))) * (4*exp(1)*nuf + 4*exp(1+rho*t)*nuf + exp(2/(1+exp(rho*t)) + rho*t)*(rho - 4*nuf) - exp(2/(1+exp(rho*t)))*(4*nuf+rho))) / (4*(1+exp(rho*t))*n*rho) );
  }

  if(output == 4)
  {
    return(sigma*sigma)
  }
  if(output == 5)
  {
    return( log((mu*(0.5-1))/(0.5*(mu-1))) )
  }
  
  dist.mean = log(mu/(1-mu)) + (2*mu-1)*sigma*sigma / (2*(mu-1)*(mu-1)*mu*mu)
  dist.var = sigma*sigma / ((mu-1)*(mu-1)*mu*mu)
  dist.sd = sqrt(dist.var)
  
  # what to output? negative log likelihood, mean, or sd?
  if(output == 1)
  {
    # log likelihood of our data
    logl = sum(dnorm(Y, dist.mean, dist.sd, TRUE))
    return(-logl) 
  }
  else if(output == 2)
  {
    return(dist.mean)
  }
  else if(output == 3)
  {
    return(dist.sd)
  }
}

# get data
df = read.table("Data/le-data.txt")

# find max lik version of models with and without each parameter
ics = c(log(0.75), log(0.01), -0.0005, log(1000))
opt.model.0 = optim(ics, model.lik, t=df[,2], Y=df[,4], model=0)
opt.model.1 = optim(ics, model.lik, t=df[,2], Y=df[,4], model=1)
opt.model.2 = optim(ics, model.lik, t=df[,2], Y=df[,4], model=2)
opt.model.3 = optim(ics, model.lik, t=df[,2], Y=df[,4], model=3)
opt.model.4 = optim(ics, model.lik, t=df[,2], Y=df[,4], model=4)
opt.model.5 = optim(ics, model.lik, t=df[,2], Y=df[,4], model=5)

# get AICs
aic.0 = 2*2 - 2*(-opt.model.0$value)
aic.1 = 2*1 - 2*(-opt.model.1$value)
aic.2 = 2*1 - 2*(-opt.model.2$value)
aic.3 = 2*4 - 2*(-opt.model.3$value)
aic.4 = 2*3 - 2*(-opt.model.4$value)
aic.5 = 2*3 - 2*(-opt.model.5$value)

# get best parameters
theta.0 = opt.model.0$par
theta.1 = opt.model.1$par
theta.2 = opt.model.2$par
theta.3 = opt.model.3$par
theta.4 = opt.model.4$par
theta.5 = opt.model.5$par

vprime.tab = read.table("Data/le-age-vprime.csv")

preds.0 = model.lik(theta.0, vprime.tab[,1], 0, model=0, output=4)
preds.1 = model.lik(theta.1, vprime.tab[,1], 0, model=1, output=4)
preds.2 = model.lik(theta.2, vprime.tab[,1], 0, model=2, output=4)
preds.3 = model.lik(theta.3, vprime.tab[,1], 0, model=3, output=4)
preds.4 = model.lik(theta.4, vprime.tab[,1], 0, model=4, output=4)
preds.5 = model.lik(theta.5, vprime.tab[,1], 0, model=5, output=4)

t.seq = seq(0, 300)
preds.0.line = model.lik(theta.0, t.seq, 0, model=0, output=4)
preds.1.line = model.lik(theta.1, t.seq, 0, model=1, output=4)
preds.2.line = model.lik(theta.2, t.seq, 0, model=2, output=4)
preds.3.line = model.lik(theta.3, t.seq, 0, model=3, output=4)
preds.4.line = model.lik(theta.4, t.seq, 0, model=4, output=4)
preds.5.line = model.lik(theta.5, t.seq, 0, model=5, output=4)

r2.0 = 1 - sum((preds.0-vprime.tab[,2])**2) / sum((vprime.tab[,2]-mean(vprime.tab[,2]))**2)
r2.1 = 1 - sum((preds.1-vprime.tab[,2])**2) / sum((vprime.tab[,2]-mean(vprime.tab[,2]))**2)
r2.2 = 1 - sum((preds.2-vprime.tab[,2])**2) / sum((vprime.tab[,2]-mean(vprime.tab[,2]))**2)
r2.3 = 1 - sum((preds.3-vprime.tab[,2])**2) / sum((vprime.tab[,2]-mean(vprime.tab[,2]))**2)
r2.4 = 1 - sum((preds.4-vprime.tab[,2])**2) / sum((vprime.tab[,2]-mean(vprime.tab[,2]))**2)
r2.5 = 1 - sum((preds.5-vprime.tab[,2])**2) / sum((vprime.tab[,2]-mean(vprime.tab[,2]))**2)

r2.log.0 = 1 - sum((log(preds.0)-log(vprime.tab[,2]))**2) / sum((log(vprime.tab[,2])-mean(log(vprime.tab[,2])))**2)
r2.log.1 = 1 - sum((log(preds.1)-log(vprime.tab[,2]))**2) / sum((log(vprime.tab[,2])-mean(log(vprime.tab[,2])))**2)
r2.log.2 = 1 - sum((log(preds.2)-log(vprime.tab[,2]))**2) / sum((log(vprime.tab[,2])-mean(log(vprime.tab[,2])))**2)
r2.log.3 = 1 - sum((log(preds.3)-log(vprime.tab[,2]))**2) / sum((log(vprime.tab[,2])-mean(log(vprime.tab[,2])))**2)
r2.log.4 = 1 - sum((log(preds.4)-log(vprime.tab[,2]))**2) / sum((log(vprime.tab[,2])-mean(log(vprime.tab[,2])))**2)
r2.log.5 = 1 - sum((log(preds.5)-log(vprime.tab[,2]))**2) / sum((log(vprime.tab[,2])-mean(log(vprime.tab[,2])))**2)

plot(vprime.tab[,1], log(vprime.tab[,2]))
points(t.seq, log(preds.3.line), type="l")
points(t.seq, log(preds.4.line), type="l")
points(t.seq, log(preds.5.line), type="l")
#points(t.seq, log(preds.3.line), type="l")

mean.tab = read.table("Data/le-age-mean.csv")
preds.mean.0 = model.lik(theta.0, mean.tab[,1], 0, model=0, output=5)
preds.mean.3 = model.lik(theta.3, mean.tab[,1], 0, model=3, output=5)
r2.mean.0 = 1 - sum((preds.mean.0-mean.tab[,2])**2) / sum((mean.tab[,2]-mean(mean.tab[,2]))**2)
r2.mean.3 = 1 - sum((preds.mean.3-mean.tab[,2])**2) / sum((mean.tab[,2]-mean(mean.tab[,2]))**2)

# construct traces of moments with best parameterisations
#for(t in 0:400) {
#  mean.0 = model.lik(theta.0, t, 0, model = 0, output = 2)
#  sd.0 = model.lik(theta.0, t, 0, model = 0, output = 3)
#  mean.1 = model.lik(theta.1, t, 0, model = 1, output = 2)
#  sd.1 = model.lik(theta.1, t, 0, model = 1, output = 3)
#  if(t == 0) {
#    moments.table = c(t, mean.0, sd.0, mean.1, sd.1)
#  }
#  else {
#    moments.table = rbind(moments.table, c(t, mean.0, sd.0, mean.1, sd.1))
#  }
#}
#write.table(moments.table, file="models-selection.txt", row.names=F, col.names=F)

