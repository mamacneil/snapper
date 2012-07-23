snapper
=======

PyMC implementation of Smith et al. (2012) linked occupancy-abundance model

Created by M. Aaron MacNeil on 20/07/12.

"Smith, Adam Nicholas Howard, Marti Jane Anderson, and Russell Brian Millar. In press. Incorporating the intraspecific occupancy-abundance relationship into zero-inflated models. Ecology. http://dx.doi.org/10.1890/12-0460.1
Incorporating the intraspecific occupancy-abundance relationship into zero-inflated models

Adam Nicholas Howard Smith1,*, Marti Jane Anderson2, and Russell Brian Millar3
1Massey University, Albany Campus, New Zealand Institute for Advanced Study

2Massey University, New Zealand Institute for Advanced Study

3University of Auckland, Department of Statistics

Zero-inflated versions of standard distributions for count data are often required in order to account for excess zeros when modeling the abundance of organisms. Such distributions typically have as parameters λ, the mean of the count distribution, and π, the probability of an excess zero. Implementations of zero-inflated models in ecology typically model λ using a set of predictor variables, and π is fit either as a constant or with its own separate model. Neither of these approaches makes use of any relationship that might exist between π and λ. However, for many species, the rate of occupancy is closely and positively related to its average abundance. Here, this relationship was incorporated into the model by functionally linking π to λ, and demonstrated in a study of snapper (Pagrus auratus) in and around a marine reserve. This approach has several potential practical advantages, including better computational performance and more straightforward model interpretation. It is concluded that, where appropriate, directly linking π to λ can produce more ecologically accurate and parsimonious statistical models of species abundance data."


The Smith et al. code OpenBUGS code is below, from their supplemental information. This PyMC version is not identical to their approach (but close enough).




model{
  for(i in 1:N) {
    iz[i] <- step(sna[i] - 1)
    NBProbZero[i] <- pow(1 + lambda[i] / delta, -delta)
    nllZero[i] <- -log(p[i] + (1 - p[i])*NBProbZero[i])
    llNB[i] <- (loggam(delta + sna[i]) - loggam(delta) - loggam(sna[i] + 1) + sna[i] * log(lambda[i]) - sna[i] * log(delta) - (delta + sna[i]) * log(1 + lambda[i] / delta))
    nllNonZero[i] <- -(log(1 - p[i]) + llNB[i])	
    nll[i] <- iz[i] * nllNonZero[i] + (1 - iz[i]) * nllZero[i]	
    zero[i] <- 0
    zero[i] ~ dpois(nll[i]) 
  
    logit(p[i]) <- gamma0 + gamma1*eta[i]
    log(lambda[i]) <- eta[i]
    eta[i] <- llam.area[area[i]] + beta.season*season[i] + beta.year[year[i]]
  }
    
  # priors  
  delta ~ dgamma(0.0001,0.0001)
  gamma0 ~ dnorm(0.0, 0.01)
  gamma1 ~ dnorm(0.0, 0.01)

  for(j in 1:6) { llam.area[j] <- llam.reserve[reserve.in.area[j]] + beta.area[j]
                  beta.area[j] ~ dnorm(0.0, tau.area)  }

  llam.reserve[1] <- beta0 - beta.reserve
  llam.reserve[2] <- beta0 + beta.reserve
	
  beta0 ~ dnorm(0.0, 0.01)
  beta.reserve ~ dnorm(0.0, 0.01)
  beta.season ~ dnorm(0.0, 0.01)
  	
  for(m in 1:9) { beta.year[m] ~ dnorm(0.0,tau.year) }
	
  # hyperpriors
  tau.year <- pow(sd.year, -2)
  sd.year ~ dunif(0,5)
  nlh.year <- log(1 + sd.year*sd.year)
  z.year <- 0
  z.year ~ dpois(nlh.year)

  tau.area <- pow(sd.area, -2)
  sd.area ~ dunif(0,5)
  nlh.area <- log(1 + sd.area*sd.area)
  z.area <- 0
  z.area ~ dpois(nlh.area)

  ## CALCULATIONS

  # variance components for fixed effects
  sd.reserve <- sd(llam.reserve[])
  sd.season <- sqrt( 2* pow(beta.season, 2) )

  # Means for no-reserve and reserve
  p.nores <- 1 / (1 + exp( -( gamma0 + gamma1*llam.reserve[1] )) )
  p.res <-  1 / (1 + exp( -( gamma0 + gamma1*llam.reserve[2] )) )
  mean.nores <- (1 - p.nores) * exp(llam.reserve[1])
  mean.res <- (1 - p.res) * exp(llam.reserve[2])
  res.effect <- mean.res / mean.nores
}



