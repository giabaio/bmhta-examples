# Loads relevant packages
library(tidyverse)
library(R2jags)

# Loads the new module
load.module("pcprob",path=".")

# Sets data
y=4
n=6
# lambda>0 controls how strictly the PC concentrates around the base-value theta0
#   higher lambda implies larger mass around theta0, while smaller lambda 
#   indicates that the mass is spread out more widely
lambda=1
theta0=0.5

# 1. Uses the JAGS module
model=function(){
  theta~dpcprob(lambda,theta0)
  y~dbin(theta,n)
}
m1=jags(
  data=list(y=y,n=n,lambda=lambda,theta0=theta0),inits=NULL,
  parameters.to.save="theta",
  model.file=model,n.chains=2,n.iter=100000,n.thin=4,DIC=FALSE
)

# 2. Uses zero-tricks
model_pc_bern=function(){
  # "Pseudo-observation", w=0
  w ~ dpois(phi)
  # Defines the likelihood from the PC prior (guards against d=0 at theta=theta0)
  d <- pow(
    2*theta*log(theta/theta0)+2*(1-theta)*log((1-theta)/(1-theta0)), 0.5
  ) + 1e-10  
  deriv <- abs((log(theta)-log(1-theta)-log(theta0)+log(1-theta0))/d)
  # Models phi as a function of the likelihood for theta
  phi <- -log(lambda*exp(-lambda*d)*deriv) + C
  # "Non-informative" prior for theta, to combine with the likelihood
  # contribution from phi and effectively imply that theta ~ PC prior
  theta ~ dbeta(1,1)
  # Possible data - if not observed, sampling from the prior only
  y ~ dbin(theta,n)
}
# Runs the model in JAGS with no data for y to check what the prior looks like
m2=jags(
  data=list(y=y,n=n,lambda=lambda,theta0=theta0,w=0,C=10000),
  parameters.to.save=c("theta"),model.file=model_pc_bern,n.chains=2,
  n.thin=4,n.iter=100000,DIC=FALSE
)

# Shows the results
print(m1)
print(m2)

# Compares the densities
bmhe::posteriorplot(m1)$data |> mutate(model="JAGS module") |> 
  bind_rows(bmhe::posteriorplot(m2)$data |> mutate(model="Zero-trick")) |> 
  ggplot2::ggplot(aes(x=value,col=model)) + geom_density(key_glyph="path") + 
  theme_bw() + theme(legend.position="bottom") + xlab("theta") + ylab("Density") + 
  labs(
    title=paste0(
      "PC Prior(lambda=",lambda,", theta_0=",theta0,"). Observed data: y=",y,", n=",n
    )
  )

# Compares the 95% estimates for the model parameter
bmhe::coefplot(m1)$data |> mutate(model="JAGS module") |> 
  bind_rows(bmhe::coefplot(m2)$data |> mutate(model="Zero-trick")) |> 
  bmhe::coefplot(xintercept=NULL) + aes(color=model) + 
  theme(legend.position="bottom")

# Tail area probabilities
sum(m1$BUGSoutput$sims.list$theta>0.75)/m1$BUGSoutput$n.sims
sum(m2$BUGSoutput$sims.list$theta>0.75)/m2$BUGSoutput$n.sims

# Can compare analysis of observed data given different priors
model_uniform=function(){
  # Uniform prior
  theta ~ dbeta(1,1)
  y ~ dbin(theta,n)
}
model_jeffreys=function(){
  # Jeffreys' prior
  theta ~ dbeta(.5,.5)
  y ~ dbin(theta,n)
}
m_uniform=jags(
  data=list(y=y,n=n),inits=NULL,
  parameters.to.save="theta",
  model.file=model_uniform,n.chains=2,n.iter=100000,n.thin=4,DIC=FALSE
)
m_jeffreys=jags(
  data=list(y=y,n=n),inits=NULL,
  parameters.to.save="theta",
  model.file=model_jeffreys,n.chains=2,n.iter=100000,n.thin=4,DIC=FALSE
)

# Visualise the output
bmhe::posteriorplot(m1)$data |> mutate(model="PC prior") |> 
  bind_rows(bmhe::posteriorplot(m_uniform)$data |> mutate(model="Uniform prior")) |> 
  bind_rows(bmhe::posteriorplot(m_jeffreys)$data |> mutate(model="Jeffreys' prior")) |> 
  ggplot2::ggplot(aes(x=value,col=model)) + geom_density(key_glyph="path") + 
  theme_bw() + theme(legend.position="bottom") + xlab("theta") + ylab("Density") 

# Compares the 95% estimates for the model parameter
bmhe::coefplot(m1)$data |> mutate(model="PC prior") |> 
  bind_rows(bmhe::coefplot(m_uniform)$data |> mutate(model="Uniform prior")) |> 
  bind_rows(bmhe::coefplot(m_jeffreys)$data |> mutate(model="Jeffreys' prior")) |> 
  bmhe::coefplot(xintercept=NULL) + aes(color=model) + 
  theme(legend.position="bottom")
