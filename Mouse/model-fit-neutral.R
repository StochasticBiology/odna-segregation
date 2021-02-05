model.lik <- function(theta, t, Y, model = 0, output=1) {

  # grab the model variables from the parameter vector
  # we use logged values to ensure they all remain non-negative
  nuf = exp(theta[1])
  v0 = exp(theta[2])

  mu = 1/2+(t-t)
  if(model == 0) {
    sigma = sqrt(v0 + nuf*2*t)
  }
  else if(model == 1) {
    sigma = sqrt(nuf*2*t)
  }
  else if(model == 2) {
    sigma = sqrt(v0)
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
df = read.table("Data/hb-data.txt")

# find max lik version of models with and without each parameter
opt.model.0 = optim(log(c(1, 10)), model.lik, t=df[,2], Y=df[,4], model=0)
opt.model.1 = optim(log(c(1, 10)), model.lik, t=df[,2], Y=df[,4], model=1)
opt.model.2 = optim(log(c(1, 10)), model.lik, t=df[,2], Y=df[,4], model=2)

# get AICs
aic.0 = 2*2 - 2*(-opt.model.0$value)
aic.1 = 2*1 - 2*(-opt.model.1$value)
aic.2 = 2*1 - 2*(-opt.model.2$value)

# get best parameters
theta.0 = opt.model.0$par
theta.1 = opt.model.1$par
theta.2 = opt.model.2$par

# construct traces of moments with best parameterisations
for(t in 0:400) {
  mean.0 = model.lik(theta.0, t, 0, model = 0, output = 2)
  sd.0 = model.lik(theta.0, t, 0, model = 0, output = 3)
  mean.1 = model.lik(theta.1, t, 0, model = 1, output = 2)
  sd.1 = model.lik(theta.1, t, 0, model = 1, output = 3)
  mean.2 = model.lik(theta.2, t, 0, model = 2, output = 2)
  sd.2 = model.lik(theta.2, t, 0, model = 2, output = 3)
  if(t == 0) {
    moments.table = c(t, mean.0, sd.0, mean.1, sd.1, mean.2, sd.2)
  }
  else {
    moments.table = rbind(moments.table, c(t, mean.0, sd.0, mean.1, sd.1, mean.2, sd.2))
  }
}
write.table(moments.table, file="models-neutral.txt", row.names=F, col.names=F)

# bootstrap to get distributions on parameters
nuf.vals = v0.vals = NULL
for(boot in 1:500) {
  r = sample.int(length(df[,1]), length(df[,1]), replace=T)
  t.boot = df[r,2]
  Y.boot = df[r,4]
  opt.model.0.boot = optim(log(c(1, 10)), model.lik, t=t.boot, Y=Y.boot, model=0)  
  nuf.vals = c(nuf.vals, exp(opt.model.0.boot$par[1]))
  v0.vals = c(v0.vals, exp(opt.model.0.boot$par[2]))
}

scatter.points = cbind(nuf.vals, v0.vals)
write.table(scatter.points, "nuf-v0-scatter.txt", row.names=F, col.names=F)

# make histograms of these distributions
nuf.breaks = seq(0,0.0002,0.00001)
v0.breaks = seq(0,0.021,0.001)
nuf.bins = as.numeric(paste(summary(cut(nuf.vals, breaks=nuf.breaks))))
v0.bins = as.numeric(paste(summary(cut(v0.vals, breaks=v0.breaks))))
nuf.hist = cbind(nuf.breaks, nuf.bins/sum(nuf.bins))
v0.hist = cbind(v0.breaks, v0.bins/sum(v0.bins))
write.table(nuf.hist, "nuf-hist.txt", row.names=F, col.names=F)
write.table(v0.hist, "v0-hist.txt", row.names=F, col.names=F)
