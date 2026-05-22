## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| eval: true
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


## ----pvalue--------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: true
pt(q=t,df=n-1,lower.tail=FALSE)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ttestFish2
#| fig-cap: 'The sampling distribution for the statistic $T$ under the null hypothesis, $p(T \mid H_0)$. The shaded area indicates the $p-$value, $\Pr(T>t=2.45\mid H_0)=0.0077$'
#| dev: "tikz"
#| echo: false
ggplot(NULL, aes(c(-4,4))) +
  geom_area(stat="function",fun=dt,args=list(df=n-1),fill="white",xlim=c(-4,t),col="black",alpha=0.15) +
  geom_area(stat="function",fun=dt,args=list(df=n-1),fill="grey",xlim=c(t,4),col="black") + 
  theme_bw() + xlab("$T$") + ylab("Sampling distribution for $T$") + 
  annotate("text",x=t,y=-.005,label=paste("$t=",format(t,digits=3,nsmall=2),"$"))

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| dev: "tikz"
#| label: fig-intervals-experts
#| fig-cap: "Even with very different starting points, the two experts become increasingly aligned in their (posterior) judgement by effect of increasingly large and definitive evidence from the data"
#| fig-height: 4.5
#| fig-width: 6

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

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-covid-testing-2
#| fig-pos: "h"
#| fig-cap: "The relationship between the 'posterior' probability of disease given the negative result and the 'prior' probability of the disease given current knowledge, considering a range of possible values in $[0;1]$"
#| dev: "tikz"
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


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Simulates a sample of S values from y ~ Normal(0,1)
S=100000
y=rnorm(S,0,1)

# Computes the 2.5%  and 97.5% quantile; these are the points y_{0.025} and 
# y_{0.975} such that Pr(Y<=y_{0.025})=0.025 and Pr(Y<=y_{0.975})=0.975
quantile(y,c(0.025,0.975))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| dev: "tikz"
#| label: fig-monte-carlo-approx
#| fig-cap: "Monte Carlo approximation of the 2.5% and 97.5% percentiles of a Normal(0,1) distribution. Increasing number of MC samples improves the accuracy and, effectively, reduces the approximation error to virtually nothing. The dashed horizontal lines indicate the analytic values for the 2.5% and the 97.5% quantiles of the standard Normal distribution"
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

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-prior-mu
#| dev: "tikz"
#| out-width: "100%"
#| fig-cap: "A prior for $\\mu\\sim\\mbox{log-Normal}(5.2,0.2)$"
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


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-prior-sigma
#| out-width: "100%"
#| dev: "tikz"
#| fig-cap: "A prior for $\\sigma\\sim\\mbox{Exponential}(0.35)$"
tibble(sigma=seq(0,15,.1)) |> mutate(p2=dexp(seq(0,15,.1),.35)) |> 
  ggplot(aes(sigma,p2)) + 
  geom_line() + xlab("$\\sigma$") + ylab("")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| out-width: "100%"
#| label: fig-prior-eta
#| dev: "tikz"
#| fig-cap: "The implied prior for $\\nu=\\mu\\gamma$"
mu=rlnorm(10000,5.2,.2)
sigma=rexp(10000,.35)
gamma=sqrt(mu/sigma^2)
nu=mu*gamma
tibble(nu) |> ggplot(aes(nu)) + geom_density() + 
  xlab("$\\nu$") + ylab("") + xlim(0,30000)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| out-width: "100%"
#| label: fig-prior-lambda
#| dev: "tikz"
#| fig-cap: "The implied prior for $\\displaystyle\\gamma=\\sqrt{\\frac{\\mu}{\\sigma^2}}$"
tibble(gamma) |> ggplot(aes(gamma)) + geom_density() + 
  xlab("$\\gamma$") + ylab("")+ xlim(c(0,200))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: true
# Defines the number of simulations
nsim = 10000  

# Then defines the prior distribution for theta
alpha = 9.2
beta = 13.8
theta = rbeta(n=nsim, shape1=alpha, shape2=beta)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: true
# Produces summary statistics for the prior
bmhe::stats(theta)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| dev: "tikz"
#| label: fig-hist-drug-prior
#| fig-cap: "A histogram describing the Monte Carlo simulations from a Beta(9.2, 13.8) distribution. As is possible to see, the distribution is roughly speaking centered around 0.4 and most values are included in the interval [0.2; 0.6], as required"
#| fig-pos: "h"
tibble(theta) |> ggplot(aes(theta)) + geom_histogram(col="black",fill="grey") +  
  xlab("$\\theta$") 


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: true
sum(theta>0.5)/length(theta)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-info-distr
#| fig-cap: "The prior *information* indicates that most of the mass for the parameter $\\theta$ should be contained between 0.2 and 0.6, which can be encoded into two different prior *distributions*: the curve shows a $\\mbox{Beta}(9,2, 13.8)$, while the histogram shows a $\\mbox{Normal}(-0.4, 0.413)$ on the logit scale. For all intents and purposes, the two are effectively identical"
#| dev: "tikz"
#| fig-pos: "h"
#| fig-height: 3.5
#| fig-width: 5.5
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


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
# Need to re-run this chunk as it can't be ported from a different file...
# May need to create a .qmd or .R file that has these instructions to run 
# once every new chapter if needed throughout...
# Defines the number of simulations
nsim = 10000  

# Then defines the prior distribution for theta
alpha = 9.2
beta = 13.8
theta = rbeta(n=nsim, shape1=alpha, shape2=beta)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
dlogitnorm=function(x,mu,sigma) {
  dnorm(log(x/(1-x)),mu,sigma)*abs(1/(x*(1-x)))
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-beta-logitnormal
#| fig-cap: "Comparison on the Beta and the logit-Normal distributions for $\\theta$"
#| dev: "tikz"
#| fig-pos: "H"
#| fig-height: 3.5
#| fig-width: 5.5
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


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
dlogitnorm(0.4,-0.405,0.413)
dbeta(0.4,9.2,13.8)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: true
# Simulates from the (prior) predictive distribution of 'y'
y=rbinom(n=nsim,size=20,prob=theta)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: true
P.crit=(y>=15)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-predictive-drug
#| fig-cap: 'The (prior) predictive distribution $p(y)$. The bars in lighter colour indicate the "tail area probability" that the observed results would exceed the threshold set at at least 15 successes out of 20 individuals'
#| dev: "tikz"
tibble(y=y) |> mutate(P.crit=y>=15) |> ggplot(aes(y,fill=P.crit)) + 
  geom_bar(stat="count",col="black") + xlab("$y$") + 
  theme(
    legend.position="none"
  ) + scale_fill_manual(values = c("#1F77B4","#FF7F0E"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: true
# Binds the vectors 'P.crit', 'y' and 'theta' into columns of a matrix
sims = cbind(P.crit,y,theta)
# Runs bmhe::stats to summarise the simulations
bmhe::stats(sims)

