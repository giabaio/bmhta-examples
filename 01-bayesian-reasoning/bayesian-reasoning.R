#' This script contains all the code used in Chapter 1. Graphs may look 
#' different than the published version, especially if they contain some 
#' mathematical notation, processed using LaTeX. If the graph is rendered 
#' using the `tikzDevice` package, then the mathematical symbols (e.g. $x$) are
#' processed correctly through LaTeX and display as in the published version.
#' In pure R, they are considered as verbatim text string.

## ----Load relevant packages----------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

## ----Set up significance test--------------------------------------------------------------------------------------------------------------------------------------
set.seed(120)
n0=80; n1=78; n=n0+n1
mu0=25; mu1=32
sigma0=17; sigma1=20
y0=rnorm(n0,mu0,sigma0)
y1=rnorm(n1,mu1,sigma1)
y0bar=mean(y0)
y1bar=mean(y1)
s2_0=sum((y0-y0bar)^2)/(n0-1) # same as var(y0)
s2_1=sum((y1-y1bar)^2)/(n1-1) # same as var(y1)
d=y1bar-y0bar
s2_D=(s2_0/n0)+(s2_1/n1)
delta=0
t=(d-delta)/sqrt(s2_D)

pt(q=t,df=n-1,lower.tail=FALSE)

ggplot(NULL, aes(c(-4,4))) +
  geom_area(stat="function",fun=dt,args=list(df=n-1),fill="white",xlim=c(-4,t),col="black",alpha=0.15) +
  geom_area(stat="function",fun=dt,args=list(df=n-1),fill="grey",xlim=c(t,4),col="black") + 
  theme_bw() + xlab("$T$") + ylab("Sampling distribution for $T$") + 
  annotate("text",x=t,y=-.005,label=paste("$t=",format(t,digits=3,nsmall=2),"$"))

## ----Experts-------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1234)
n=10
y=sum(rbinom(n,1,0.93))
a0=3; b0=12; c0=12; d0=3
a1=a0+y; b1=b0+n-y
c1=c0+y; d1=d0+n-y
n1=200
y1=sum(rbinom(n1,1,0.93))
a2=a0+y1; b2=b0+n1-y1
c2=c0+y1; d2=d0+n1-y1

ints=tibble(
  lower=c(
    rbeta(100000,a0,b0) |> quantile(0.025),
    rbeta(100000,c0,d0) |> quantile(0.025),
    rbeta(100000,a1,b1) |> quantile(0.025),
    rbeta(100000,c1,d1) |> quantile(0.025),
    rbeta(100000,a2,b2) |> quantile(0.025),
    rbeta(100000,c2,d2) |> quantile(0.025)
  ),
  upper=c(
    rbeta(100000,a0,b0) |> quantile(0.975),
    rbeta(100000,c0,d0) |> quantile(0.975),
    rbeta(100000,a1,b1) |> quantile(0.975),
    rbeta(100000,c1,d1) |> quantile(0.975),
    rbeta(100000,a2,b2) |> quantile(0.975),
    rbeta(100000,c2,d2) |> quantile(0.975)
  ),
  expert=rep(c("$E_1$","$E_2$"),3),
  time=factor(
    rep(c("Prior",paste0("Posterior with\n $y_1=",y,",n_1=",n,"$"),paste0("Posterior with\n $y_2=",y1,",n_2=",n1,"$")),each=2),
    level=c(paste0("Posterior with\n $y_2=",y1,",n_2=",n1,"$"),paste0("Posterior with\n $y_1=",y,",n_1=",n,"$"),"Prior")
  )
) 

ints |> 
  ggplot() + xlim(0,1) + xlab("Interval estimate for $\\theta$") + ylab("") +
  geom_segment(aes(x=lower,xend=upper,y=time,yend=time,col=expert),linewidth=1.1) +
  geom_segment(
    data=tibble(prior=ints$lower[1],posterior=ints$lower[3]),aes(x=prior,xend=posterior,y=3,yend=2),
    arrow = arrow(length=unit(0.20,"cm"),type="closed"),col="gray70",size=.2
  ) +
  geom_segment(
    data=tibble(prior=ints$lower[2],posterior=ints$lower[4]),aes(x=prior,xend=posterior,y=3,yend=2),
    arrow = arrow(length=unit(0.20,"cm"),type="closed"),col="gray70",size=.2
  ) +
  geom_segment(
    data=tibble(prior=ints$upper[1],posterior=ints$upper[3]),aes(x=prior,xend=posterior,y=3,yend=2),
    arrow = arrow(length=unit(0.20,"cm"),type="closed"),col="gray70",size=.2
  ) +
  geom_segment(
    data=tibble(prior=ints$upper[2],posterior=ints$upper[4]),aes(x=prior,xend=posterior,y=3,yend=2),
    arrow = arrow(length=unit(0.20,"cm"),type="closed"),col="gray70",size=.2
  ) +
  geom_segment(
    data=tibble(prior=ints$lower[3],posterior=ints$lower[5]),aes(x=prior,xend=posterior,y=2,yend=1),
    arrow = arrow(length=unit(0.20,"cm"),type="closed"),col="gray70",size=.2
  ) +
  geom_segment(
    data=tibble(prior=ints$lower[4],posterior=ints$lower[6]),aes(x=prior,xend=posterior,y=2,yend=1),
    arrow = arrow(length=unit(0.20,"cm"),type="closed"),col="gray70",size=.2
  ) +
  geom_segment(
    data=tibble(prior=ints$upper[3],posterior=ints$upper[5]),aes(x=prior,xend=posterior,y=2,yend=1),
    arrow = arrow(length=unit(0.20,"cm"),type="closed"),col="gray70",size=.2
  ) +
  geom_segment(
    data=tibble(prior=ints$upper[4],posterior=ints$upper[6]),aes(x=prior,xend=posterior,y=2,yend=1),
    arrow = arrow(length=unit(0.20,"cm"),type="closed"),col="gray70",size=.2
  ) + theme_classic() +
  theme(
    legend.title=element_blank(),
    legend.position=c(0.15, 0.15),
    legend.background=element_rect(fill='transparent'),
    axis.line.y=element_blank(),axis.ticks.y=element_blank()
  ) + scale_color_manual(values = c("#E69F00","#0072B2"))


## ----Covid testing-------------------------------------------------------------------------------------------------------------------------------------------------
theta=seq(0,1,.01)
data=tibble(prior=theta, post=theta*0.04/(theta*0.04 + (1-theta)*.95))
data %>% ggplot() + geom_line(aes(prior,post)) + 
  labs(x="Prior probability of disease, $\\theta\\mid\\mathcal{B}$", y="Posterior given -ve") +
  annotate("text",.1,data %>% dplyr::filter(prior==0.1) %>% pull(post),label=format(data %>% dplyr::filter(prior==0.1) %>% pull(post),digits=4),vjust=-1.2) +
  annotate("text",.4,data %>% dplyr::filter(prior==0.4) %>% pull(post),label=format(data %>% dplyr::filter(prior==0.4) %>% pull(post),digits=4),vjust=-1.2) +
  annotate("text",.8,data %>% dplyr::filter(prior==0.8) %>% pull(post),label=format(data %>% dplyr::filter(prior==0.8) %>% pull(post),digits=4),vjust=-1.2) +
  geom_segment(aes(x=.1,y=-Inf,xend=.1,yend=data %>% dplyr::filter(prior==0.1) %>% pull(post)), linetype="dashed") +
  geom_segment(aes(x=.4,y=-Inf,xend=.4,yend=data %>% dplyr::filter(prior==0.4) %>% pull(post)), linetype="dashed") +
  geom_segment(aes(x=.8,y=-Inf,xend=.8,yend=data %>% dplyr::filter(prior==0.8) %>% pull(post)), linetype="dashed") 


## ----Monte Carlo approximation-------------------------------------------------------------------------------------------------------------------------------------
# Simulates a sample of S values from y ~ Normal(0,1)
S=100000
y=rnorm(S,0,1)

# Computes the 2.5%  and 97.5% quantile; these are the points y_{0.025} and 
# y_{0.975} such that Pr(Y<=y_{0.025})=0.025 and Pr(Y<=y_{0.975})=0.975
quantile(y,c(0.025,0.975))

S=c(10,15,20,30,50,75,100,200,300,500,750,1000,2000,5000,7500,10000,50000,100000,500000,1000000)
y=lapply(S,function(S) rnorm(S,0,1))
low=lapply(1:length(S),function(i) quantile(y[[i]],0.025)) |> unlist()
upp=lapply(1:length(S),function(i) quantile(y[[i]],0.975)) |> unlist()
tibble(S=S,quantile=low,type="0.025 quantile") |> 
  bind_rows(
    tibble(S=S,quantile=upp,type="0.975 quantile")
  ) |> ggplot(aes(S,quantile,col=type)) + geom_point() + 
  geom_line() +
  geom_hline(yintercept=qnorm(0.025),linetype=2) + 
  geom_hline(yintercept=qnorm(0.975),linetype=2) + 
  xlab("Number of Monte Carlo samples, $S$") + ylab("Monte Carlo estimates") +
  theme(
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.position = c(.75,.5)
  ) + scale_x_continuous(
    trans='log',labels=scales::comma,
    breaks=c(10,100,1000,10000,100000,1000000)
  ) + scale_color_manual(values = c("#E69F00","#0072B2"))


## ----Forward sampling----------------------------------------------------------------------------------------------------------------------------------------------
set.seed(3020)
# Simulate values from the prior on the *natural* parameters
mu=rlnorm(10000,5.2,.2)
sigma=rexp(10000,.35)

# Transform the prior for the *original scale* parameters
gamma=sqrt(mu/sigma^2)
nu=mu*gamma

tibble(mu=seq(0,600)) |> mutate(p1=dlnorm(seq(0,600),5.2,.2)) |> 
  ggplot(aes(mu,p1)) + 
  geom_line() + xlab("$\\mu$") + ylab("")

tibble(sigma=seq(0,15,.1)) |> mutate(p2=dexp(seq(0,15,.1),.35)) |> 
  ggplot(aes(sigma,p2)) + 
  geom_line() + xlab("$\\sigma$") + ylab("")

mu=rlnorm(10000,5.2,.2)
sigma=rexp(10000,.35)
gamma=sqrt(mu/sigma^2)
nu=mu*gamma
tibble(nu) |> ggplot(aes(nu)) + geom_density() + 
  xlab("$\\nu$") + ylab("") + xlim(0,30000)

tibble(gamma) |> ggplot(aes(gamma)) + geom_density() + 
  xlab("$\\gamma$") + ylab("")+ xlim(c(0,200))


## ----Drug----------------------------------------------------------------------------------------------------------------------------------------------------------
# Defines the number of simulations
nsim = 10000  

# Then defines the prior distribution for theta
alpha = 9.2
beta = 13.8
theta = rbeta(n=nsim, shape1=alpha, shape2=beta)

# Produces summary statistics for the prior
bmhe::stats(theta)

# Histogram
tibble(theta) |> ggplot(aes(theta)) + geom_histogram(col="black",fill="grey") +  
  xlab("$\\theta$") 

sum(theta>0.5)/length(theta)


## ----Priors---------------------------------------------------------------------------------------------------------------------------------------------------------
dat1=tibble(x=seq(0,1,.001),p1=dbeta(seq(0,1,.001),9.2,13.8))
dat2=tibble(x=rlnorm(10000,bmhe::logit(0.4),0.413))
tibble(p2=bmhe::ilogit(rnorm(10000,bmhe::logit(.4),0.413))) |>
  ggplot(aes(p2)) + geom_histogram(
    aes(
      y = stat(density),
      color="$\\mbox{log}\\left(\\frac{\\theta}{1-\\theta}\\right) \\sim \\mbox{Normal}(-0.405,0.413)$"
    ),fill="grey"
  ) +
  geom_line(
    data=dat1,aes(
      x,p1,col="$\\theta \\sim \\mbox{Beta}(9.2,13.8)$"
    ),linewidth=0.8
  ) + xlab("$\\theta$") +
  ylab("") + scale_color_manual("",values=c("black","blue")) + 
  theme(
    legend.position=c(0.75, 0.85),
    legend.background=element_rect(fill='transparent')
  )

# Defines the `dlogitnorm` function
dlogitnorm=function(x,mu,sigma) {
  dnorm(log(x/(1-x)),mu,sigma)*abs(1/(x*(1-x)))
}

# Define a dataset with x = a range of values in [0;1] and y computed as 
# their density for the logit-Normal(-0.405,0.4137) distribution
tibble(
  x=seq(0,1,.001),y=dlogitnorm(seq(0,1,.001),-.405,.4137),
  model="logit-Normal($-0.405,0.4137$)"
) |> 
# "Binds" to this dataset another one where x is the same as before and 
# y = the resulting Beta(9.2,13.8) density
  bind_rows(
    tibble(
      x=seq(0,1,.001),y=dbeta(seq(0,1,.001),9.2,13.8),
      model="Beta($9.2,13.8$)"
    )
  ) |> 
# Finally uses ggplot to make the graphical comparison
  ggplot(aes(x,y,color=model,linetype=model)) + geom_line(linewidth=0.75) +
  theme(
    legend.position="inside",legend.position.inside=c(.75,.75),
    legend.background=element_rect(fill='transparent')
  ) + 
  labs(colour="") + xlab("$\\theta$") + ylab("$p(\\theta)$") +
  scale_color_manual(values = c("#E69F00","#0072B2"), name = "") +
  scale_linetype_manual(values = c("solid","dashed"), name = "")
#  scale_color_manual(values = c("#E69F00","#0072B2"))

# Check values
dlogitnorm(0.4,-0.405,0.413)
dbeta(0.4,9.2,13.8)


## ----Predictive distribution----------------------------------------------------------------------------------------------------------------------------------------
#| echo: true
# Simulates from the (prior) predictive distribution of 'y'
y=rbinom(n=nsim,size=20,prob=theta)

P.crit=(y>=15)

tibble(y=y) |> mutate(P.crit=y>=15) |> ggplot(aes(y,fill=P.crit)) + 
  geom_bar(stat="count",col="black") + xlab("$y$") + 
  theme(
    legend.position="none"
  ) + scale_fill_manual(values = c("#1F77B4","#FF7F0E"))

# Combines the simulations
sims = cbind(P.crit,y,theta)
# Runs bmhe::stats to summarise the simulations
bmhe::stats(sims)

