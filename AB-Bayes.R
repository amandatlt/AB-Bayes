#
# Simulating an A/B test using Bayesian A/B testing techniques from Chris Stucchio, VMO
# March 2017
#
#assume values for a and b parameters for beta prior
#generate posterior probability distributions for treatment(t) and control(c)
#find joint posterior probability
#calculate expected loss
#

###---load packages---###
.libPaths("P:/User/Amanda.Teo/R/Packages")

###---generate assumptions---###

#prior
mu <- 0.1 #on average, 10% of people will convert
sigma <- 0.2 #average distance away from mu

#sample outcomes
convert_c <- 100
num_c <- 2500

convert_t <- 120
num_t <- 2500

###---visualise beta prior---###
a <- (1-mu - sigma^2*mu^3)/(sigma^2mu^2)
