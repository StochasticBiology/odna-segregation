model.lik <- function(theta, t, Y, model = 0, output=1) {

  # grab the model variables from the parameter vector
  # we use logged values to ensure they all remain non-negative
  nuf = exp(theta[1])
  v0 = exp(theta[2])
  rho = theta[3]
  n = exp(theta[4])
  
  if(model == 1) {
    mu = 1/2+(t-t)
    sigma = sqrt(v0 + nuf*2*t)
  }
  else if(model == 0) {
    mu = 1/(1 + exp(-rho*t))
    sigma = sqrt( v0 + (exp(-2/(1+exp(rho*t))) * (4*exp(1)*nuf + 4*exp(1+rho*t)*nuf + exp(2/(1+exp(rho*t)) + rho*t)*(rho - 4*nuf) - exp(2/(1+exp(rho*t)))*(4*nuf+rho))) / (4*(1+exp(rho*t)*n*rho)) );
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
opt.model.0 = optim(log(c(1, 10, 1, 1000)), model.lik, t=df[,2], Y=df[,4], model=0)
opt.model.1 = optim(log(c(1, 10, 1, 1000)), model.lik, t=df[,2], Y=df[,4], model=1)

# get AICs
aic.0 = 2*4 - 2*(-opt.model.0$value)
aic.1 = 2*2 - 2*(-opt.model.1$value)

# get best parameters
theta.0 = opt.model.0$par
theta.1 = opt.model.1$par

# construct traces of moments with best parameterisations
for(t in 0:400) {
  mean.0 = model.lik(theta.0, t, 0, model = 0, output = 2)
  sd.0 = model.lik(theta.0, t, 0, model = 0, output = 3)
  mean.1 = model.lik(theta.1, t, 0, model = 1, output = 2)
  sd.1 = model.lik(theta.1, t, 0, model = 1, output = 3)
  if(t == 0) {
    moments.table = c(t, mean.0, sd.0, mean.1, sd.1)
  }
  else {
    moments.table = rbind(moments.table, c(t, mean.0, sd.0, mean.1, sd.1))
  }
}
write.table(moments.table, file="models-selection.txt", row.names=F, col.names=F)

