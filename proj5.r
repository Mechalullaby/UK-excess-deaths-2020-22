# Likang Xu (s2295871)

# This project is about estimating the UK excess deaths 2020-22. Given the 
# observed data sets, `lt1720uk.dat`, which contains 7 columns:
# `age` - the age classes, where 0: 0-1 years
# `fpop17` and `mpop17` - the female and male populations of each age band at 
# the beginning of 2017
# `fpop20` and `mpop20` - the same but at 2020.
# `mf` and `mm` - female and male per capita death rates computed from 2017-19

# and `death1722uk.dat`, which contains:
# `deaths` - the number of deaths each week
# `week` -  the week since the start of 2017
# `d` -  the mortality rate modifier


# We will estimate the excess deaths in two ways:
# 1. Using the population model, for week j and age class i:
#             Di = 0.9885*dj*qi*Ni, 
#             Ni_star = Ni-Di, 
#             Ni_plus = Ni_star*51/52+N(i-1)_star/52, 
#             where i-1,...,101.
# Di - the number of deaths of week j
# dj - the mortality rate modifier for week j
# qi - the expected proportion of dying (qi=1-exp(-mi/52))
# mi - the annual deaths rate
# Ni - the population at the start of week j
# Ni_plus - the population at the start of week j+1
# Set N0_star = N1 at the start, given the function `pre_death`, which can 
# predict the deaths, the excess deaths from 2020 is calculated.

# 2. Using a time series model, the excess deaths in week i, xi:
#             mu1 = alpha,
#             mu(i+1) = (xi-alpha)*rho+alpha,
#             xi ~ tk(mui, tau)
#             where i=1,...,n
# tk(mui, tau) is a scaled t distribution, given the prior, tau ~ exp(1), 
# rho ~ U(0,0.9), alpha ~ N(0,tau=0.0001), k ~ U(2,100). MCMC is used the 
# draw 10000 sample from the posterior densities.


# The study ploted 6 plots to show the results, which would be explained below.



library(rjags)
library(coda)
# load the data
lt1720 <- read.table("lt1720uk.dat")
d1722 <- read.table("death1722uk.dat")


# pre_death arguments:
# Ns: the female and male populations in each 1 year age class
# ms: female and male annual death rates for each 1-year age band
# d: the mortality rate modifier vector for each week

# This function used a 156*101*2 array, `D`, to store the predicted deaths
# for each week, age class, gender. The deaths is calculated using the first 
# model. Then, returns the predicted number of deaths each week in d.

pre_death <- function(Ns, ms, d){
  # array to store the predicted total number of deaths each week each gender
  D <- array(0,dim=c(length(d),nrow(Ns),ncol(Ns)))
  # to count the total number of death for each gender
  Dt <- array(0,dim=c(length(d),ncol(Ns)))
  
  # loop of female and male
  for (z in 1:ncol(Ns)){
    N <- Ns[,z] # Ni
    m <- ms[,z] # mi
    q <- 1-exp(-m/52) # qi
    N1 <- N[1] # N0_star
    
    # loop of each week
    for (j in 1:length(d)){
      D[j,,z] <- 0.9885*d[j]*q*N # deaths
      N_star <- N-D[j,,z]
      
      # population next week
      N[1] <- N_star[1]*51/52+N1/52
      # loop of each age
      for (i in 2:length(N)){
        N[i] <- N_star[i]*51/52+N_star[i-1]/52
      }
    }
    
    # total number of death for each gender
    Dt[,z] <- rowSums(D[,,z])
  }
  
  # total number of death
  return(rowSums(Dt))
}


# total predicted deaths from the start of 2020 to the end of the data
p_t <- pre_death(lt1720[c('fpop20','mpop20')], 
                 lt1720[c('mf','mm')], d1722$d[157:length(d1722$d)])
# the observed deaths from the start of 2020 to the end of the data
o_t <- d1722$deaths[157:length(d1722$d)]

# the excess deaths from the start of 2020 to the end of the data
ex_t <- sum(o_t-p_t)


# the predicted deaths of 2020
p_20 <- pre_death(lt1720[c('fpop20','mpop20')], 
                  lt1720[c('mf','mm')], d1722$d[157:(156+156/3)])
# the observed deaths of 2020
o_20 <- d1722$deaths[157:(156+156/3)]

# the excess deaths of 2020
ex_20 <- sum(o_20-p_20)


# plot the observed deaths against week from 2020
x <- seq(1,length(p_t))
plot(x, o_t, xlab = 'Week', ylab = 'Death', ylim = c(0, (max(o_t)+2000)),
     main = paste(' The excess deaths in 2020: ', round(ex_20,2),
                  '\nThe excess deaths overall: ', round(ex_t,2)))
lines(x, p_t, col = 'blue') # plot the predicted deaths against week
legend("topright", c("observed",'predicted'), pch = c(1,NA), 
       lty = c(NA,1), col = c(1,'blue'), bty="n")


# the vector of excess deaths each week from 2020
ex_w <- o_t-p_t
# plot the cumulative excess deaths by week
plot(x, cumsum(ex_w), xlab = 'Week', ylab = 'Cumulative Excess Death',
     main = 'The cumulative excess deaths by week from 2020')


# set ex_w in weeks 51, 52, 53, 105 and 106 to NA for recording problems
backup <- ex_w[c(51,52,53,105,106)] # store the outliers for later usage
ex_w[c(51,52,53,105,106)] <- NA


# Here we construct the time series model using jags
mod <- jags.model("model.jags",data=list(x=ex_w,N=length(ex_w)))
# draw 10000 samples from the posterior densities of the mus, rho, and k
sam.coda <- coda.samples(mod, c('mu',"rho","k"), n.iter=10000)


# transform the type of the samples
sam <- as.data.frame(as.matrix(sam.coda))

# trace plot for rho
plot(sam$rho,type="l",xlab = 'Iterations', ylab=expression(rho),
     main = 'Trace of rho')

# histogram for rho
hist(sam$rho,xlab=expression(rho),main="Density of rho",probability=TRUE)


# the mean of all the parameters (mu, rho, k)
mean_all <- summary(sam.coda)$statistics[,"Mean"]
# the posterior expected value vector for mu
mean_mu <- mean_all[2:(1+149)] # since there are 149 mus


# plot every 50th sampled mu against week
plot(x, sam[50,2:(1+149)], type="l", xlab='Week', ylab=expression(mu),
     ylim=c(min(backup)-500,max(ex_w, na.rm=T)+1000), 
     main='Sampled mu', col='grey')
for (i in seq(100,10000,50)) {
  lines(x, sam[i,2:(1+149)], col='grey')
}
# plot the estimated expectation for mu
lines(x, mean_mu, col='blue')
# plot the observed excess deaths
points(x, ex_w, pch=20)
# plot the excess deaths not used for inference (outliers)
points(c(51,52,53,105,106), backup, pch=20, col='red')
legend("topright", c("every 50th mu",'mean mu','observed','outliers'), 
       pch = c(NA,NA,20,20), lty = c(1,1,NA,NA), 
       col = c('grey','blue',1,'red'), bty="n")


# plot the residuals, xi − mean_mu, against time
plot(x, ex_w-mean_mu, xlab='Week', ylab='Residuals',
     main='The residuals against time')




