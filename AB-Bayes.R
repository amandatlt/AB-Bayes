#' Simulating an A/B test using Bayesian A/B testing techniques from Chris Stucchio, VMO
#' Additional graphing techniques from ethen8181.github.io/Business-Analytics/ab_tests/frequentist_ab_test
#' March 2017
#' 
#' assume values for a and b parameters for beta prior
#' generate posterior probability distributions for treatment(t) and control(c)
#' find joint posterior probability
#' calculate expected loss


###---load packages---###
.libPaths("P:/User/Amanda.Teo/R/Packages")
library(ggplot2)

###---generate assumptions---###

#prior
mu <- 0.1 #on average, 10% of people will convert
sigma <- 0.7 #average distance away from mu

#sample outcomes
convert_c <- 20
num_c <- 100

convert_t <- 20
num_t <- 100

###---visualise beta distribution---###

a <- (1- mu - sigma^2*mu^3)/(sigma^2*mu^2)
b <- ((1 - mu)*(1-mu - sigma^2*mu^3))/(sigma^2*mu^3)

plot_beta <- function(a,b) {
  x_axis <- seq(0,1, 0.01)
  
  distribution <-data.frame(x_axis, y1 = dbeta(x_axis,a,b))
  
  plot <- ggplot(data = distribution, aes(x = x_axis)) +
          geom_line(aes(y = y1), colour = "blue", size = 1.2) +
          labs(x = "", y = "", title = sprintf("a = %.1f, b = %.1f", a, b))
        
  return(plot)        
}

###---compute posteriors for t and c given prior---###

#' posterior distribution for lamda given sample size(n) and number converted(c)
#' p(lamda|n,c) = f(x;a+c,(b+n-c)) where f is the beta distribution and a.b are parameters from the prior
#' not sure if we assume the same prior for both groups?

lamda <- seq(0,1, 0.01)

#calculate posterior of control
a_c <- a + convert_c
b_c <- b + num_c - convert_c

posterior_c <- dbeta(lamda, a_c, b_c) 

#calculate posterior of treatment in the form of  p(a < x < b)
a_t <- a + convert_t
b_t <- b + num_t - convert_t

posterior_t <- dbeta(lamda, a_t, b_t) 

###---calculate joint posterior function---###

posterior_joint <- matrix(data = NA, nrow = length(lamda), ncol = length(lamda))
colnames(posterior_joint) <- paste0("x_",lamda)
                          
for (i in 1:length(lamda)) {
  for (j in 1:length(lamda)) {
    posterior_joint[i, j] = posterior_c[i]*posterior_t[j] 
  } 
}

###---visualise joint posterior distribution---###

#place in dataframe and reshape
df <- data.frame(x_axis = seq(0,1,0.01), posterior_joint)
df <- df[rowSums(df) > 0, colSums(df) > 0]

coltoshape <- grep("0", names(df))

df<- reshape(df, idvar = "x_axis", 
             timevar = "y_axis", #name of 2nd id
             times = as.numeric(gsub("x_","",names(df)[coltoshape])), #values you want for the 2nd id
             varying = list(names(df)[coltoshape]), #in list form, columns you want to reshape
             direction = "long" )

#tidy data frame before visualisation
rownames(df) <- NULL
colnames(df)[colnames(df) != "x_axis" & colnames(df) != "y_axis"] <- "density" 
df <- df[df$density != 0,]

ggplot(data = df, aes(x = x_axis, y = y_axis)) +
  geom_point(aes(alpha = density)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  labs(x = "control conversion rate", y = "treatment conversion rate") +
  geom_abline(slope = 1, intercept = 0, colour = "blue")
  
#' shaded area is where joint probability densities of lamda_t and lamda_c are highest
#' above the blue line means treatment likely to be better

###---compute the error function---###

#' error function tells us what is the probability we make a mistake in choosing to display a particular variant
#' increased the granularity of lamda for more accurate computation

#compute here the cumlative probabilities for an interval p(a<x<b)
lamda = seq(0,1,length.out = 10001)

cumprob_c <- pbeta(lamda, a_c, b_c)
lag <- shift_by(cumprob_c, -1)
lag[is.na(lag)] <- 0
interval_c <- (cumprob_c - lag)[-1]

cumprob_t <- pbeta(lamda, a_t, b_t)
lag <- shift_by(cumprob_t, -1)
lag[is.na(lag)] <- 0
interval_t <- (cumprob_t - lag)[-1]

#calculate error function for control i.e. prob. that treatment better than control
error_treatbetter <- 0
rowtotal <- sum(interval_t)
for (i in 1:(length(lamda)-1)) {
  rowtotal <- rowtotal - interval_t[i]
  error_treatbetter <- error_treatbetter + interval_c[i]*rowtotal
}

#error function for treatment
error_controlbetter <- 0
rowtotal <- sum(interval_c)
for (i in 1:(length(lamda)-1)) {
  rowtotal <- rowtotal - interval_c[i]
  error_controlbetter <- error_controlbetter + interval_t[i]*rowtotal
}

#alternative calculation for error function of control i.e. prob that treatment is better
#NOT WORKING

f_t <- function(x) dbeta(x, a_t, b_t)
marg_probabilty_treatbetter <- function(lamda_c) {
  integrate(f_t, lower = lamda_c, upper = 1)
} 
f_c <- function(lamda_c) {
  val <- (marg_probabilty_treatbetter(lamda_c))[[1]]
  dbeta(lamda_c, a_c, b_c)*val
}
error_treatbetter_alt <- integrate(f_c, lower = 0, upper = 1)

format(error_treatbetter_alt$value, nsmall = 10)

#' note precision in probability is up to 3 decimal places only for 100,000 intervals
#' precision not great when using 10,000 intervals

###---compute the expected loss---###

#' this tells us the expected loss from choosing one variation over the other
#' first calculate loss from choosing B over A

#lamda A - lamda B (10000 x 10000 matrix)
mat_B = seq(1,1,length.out = 10001) %*% t(lamda)
mat_B = mat_B[-1,-1]
mat_A = t(mat_B)
loss_choiceB <- mat_A - mat_B

#convert loss matrix into lower triangular matrix
loss_choiceB[upper.tri(loss_choiceB, diag = TRUE)] <- 0
loss_choiceB[1:10, 1:10] #to check we did the right thing

#calculate EXPECTED loss from choosing B over A
interval_joint = interval_c %*% t(interval_t)

ptm <- proc.time()
expected_loss = 0
for (i in 1:(length(lamda)-1)) {
    expected_loss = expected_loss + t(loss_choiceB[i,]) %*% interval_joint[i,]
}
proc.time() - ptm

#currently takes about 8 seconds

#alternative method to calculate loss function
integrand <- function(x) {
  loss_func <- max(x[1] - x[2], 0)
  loss_func*dbeta(x[1], a_c, b_c)*dbeta(x[2], a_t, b_t)
}

expected_loss_alt <- adaptIntegrate(integrand, lowerLimit = c(0,0), upperLimit = c(1,1))
format(expected_loss_alt$integral, nsmall = 10)

#' expected_loss_alt gives very similar results to expected_loss when n is very large for both t and c groups
