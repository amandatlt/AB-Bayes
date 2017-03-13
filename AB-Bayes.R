#' Simulating an A/B test using Bayesian A/B testing techniques from Chris Stucchio, VMO
#' Additional graphing techniques from ethen8181.github.io/Business-Analytics/ab_tests/frequentist_ab_test
#' March 2017
#' 
#' assume values for a and b parameters for beta prior
#' generate posterior probability distributions for treatment(t) and control(c)
#' find joint posterior probability
#' calculate expected loss


###---load packages---###
#.libPaths("P:/User/Amanda.Teo/R/Packages")
library(ggplot2)

###---generate assumptions---###

#prior
mu <- 0.1 #on average, 10% of people will convert
sigma <- 0.7 #average distance away from mu

#sample outcomes
convert_c <- 10
num_c <- 100

convert_t <- 15
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

#' cannot use dbeta function because R does not compute integrals
#' \increased the granularity of lamda for more accurate computation

#compute here the cumlative probabilities for an interval p(a<x<b)
lamda = seq(0,1,length.out = 10000)

cumprob_c <- pbeta(lamda, a_c, b_c)
lag <- shift_by(cumprob_c, -1)
lag[is.na(lag)] <- 0
interval_c <- cumprob_c - lag

cumprob_t <- pbeta(lamda, a_t, b_t)
lag <- shift_by(cumprob_t, -1)
lag[is.na(lag)] <- 0
interval_t <- cumprob_t - lag

#calculate error function for treatment i.e. prob. that treatment better than control
error_treat <- 0
rowtotal <- sum(interval_t)
for (i in 2:length(lamda)) {
  rowtotal <- rowtotal - interval_t[i-1]
  error_treat <- error_treat + interval_c[i]*rowtotal
}

error_control <- 0
rowtotal <- sum(interval_c)
for (i in 2:length(lamda)) {
  rowtotal <- rowtotal - interval_c[i-1]
  error_control <- error_control + interval_t[i]*rowtotal
}

#' note precision in probability is up to 3 decimal places only

###---compute the expected loss---###

#' this tells us the expected loss from choosing one variation over the other

#loss from choosing one variation over another
loss_val <- function (lamda_a,lamda_b,choice) {
  if (choice == "A") {
    x <- max(lamda_b - lamda_a, 0)
  } else if (choice == "B") {
    x <- max(lamda_a - lamda_b,0)
  } else {
    x <- NA
  }
  return(x)
}

expected_error_A = 0
for (i in 1:length(lamda)) {
  for (j in 1:length(lamda)) {
    expected_error_A <- expected_error_A + 
                        interval_c[i]*interval_t[j]*loss_val(lamda[i], lamda[j], "A")
  }
}
 
# too slow! 