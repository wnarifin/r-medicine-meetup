# R in Medicine Meetup
# 20180809
# Wan Nor Arifin

#' # Finite mixture

# A good introduction can be found here:
# http://www.mas.ncl.ac.uk/~nmf16/teaching/mas8391/slides8303-4.pdf
# by Malcolm Farrow

#' ## About
# - Let say we have a continuous variable e.g. height that comes from
#   two groups of people. And we have no data on the grouping. Can we extract 
#   the mean height for each groups?
# - Finite mixture model (FMM) is a statistical modeling approach to do that.
# - Usually applied to continuous data.
# - In medicine, it can be used to come up with subgroups of disease.
#   e.g. severity of depression etc.
# - Latent class analysis (LCA) is the statistical analog of FMM for categorical
#   data.

#' ## Library
library(rjags)  # must also install jags in your computer before rjags

#' ## Data
# Generate data, say height in cm
# - male = N(mean = 180, SD = 7)
# - female = N(mean = 140, SD = 7)
ht = c(rnorm(150, 180, 7), rnorm(150, 140, 7))
hist(ht)
jagsdata = list(y = ht, n = length(ht))  # jagsdata must be in list
jagsdata

#' ## FMM
#' Finite mixture model,
#' $$p(\boldsymbol{y}_i|\boldsymbol\vartheta)=\sum_{k=1}^{K}\eta_kp(\boldsymbol{y}_i|\boldsymbol\theta_k)$$
#' where $\boldsymbol\theta_1, \ldots, \boldsymbol\theta_K$ = component parameters, $\boldsymbol\eta = (\eta_1, \ldots, \eta_K)$, $\boldsymbol\vartheta = (\boldsymbol\theta_1, \ldots, \boldsymbol\theta_K, \boldsymbol\eta)$ = mixture model parameters.
 
# Aim of FMM: Can we "reverse engineer" the process? i.e we obtain the male/female grouping?
#
# Suppose we assume:
# 1. two groups, k = 2
# 2. prior mean = 160, SD = 7, var = 7^2 = 49; precision, tau = 1/var = 1/49
#    mu ~ N(160, 49); in jags mu ~ N(mean, precision) ~ N(160, 1/49)
# 3. prior gamma, mean precision = 1/49, 
#    shape, a = 3; rate, b = a/mean = 3/(1/49) = 147
#    tau ~ gamma(3, 147)
# 4. prior beta, pi ~ beta(3, 3) -- gives 50:50 ratio between the groups

#' ## JAGS code
writeLines("
model {
  # data
  for (i in 1:n) {
    c[i] ~ dcat(p[1:2])
    y[i] ~ dnorm(mu[c[i]], tau[c[i]])
  }
  
  # tau
  for (j in  1:2) {
    tau[j] ~ dgamma(3, 147) # precision, tau
    sigma[j] = sqrt(1/tau[j]) # transform back to SD
  }
  
  # mu
  mu0 ~ dnorm(160, 1/50) # hierarchical mean, mu
  for (j in 1:2) {
    mu_mu0[j] ~ dnorm(mu0, 1/50) # sample mu
  }
  mu[1:2] = sort(mu_mu0)  # set order constraint, mu1 < mu2
  
  # pi
  pi ~ dbeta(3, 3)
  p[1] = pi # weight for grp 1
  p[2] = 1 - pi # weight for grp 2
}
", con = "finmix.bug")

#' ## Fit JAGS model 
jagsmodel = jags.model("finmix.bug", jagsdata, n.chains = 2)
update(jagsmodel, 1000)  # burnin = 1000
jagspars = c("mu", "sigma", "pi")
jagssamples = coda.samples(jagsmodel, jagspars, n.iter = 10000, thin = 10)
summary(jagssamples)
plot(jagssamples)
