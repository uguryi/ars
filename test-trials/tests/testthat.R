############
##  TEST  ##
############

## install.packages("testthat")
#  library(testthat)

## install.packages("rmutil")
#  library(rmutil)


## Test function
## Input:
## n: sample size
## fun: the test density function
## rfun: function that would generate n samples from the test distribution
## fun_name: name of the function
## D_left/D_right: left/right end of domain


test_fun <- function(n, fun, rfun, fun_name, D_left, D_right){
  # Generate samples using ars() and rfun()
  print(paste0("generating samples from ", fun_name,"..."))
  sample <- try(ars( g = fun, n = n, D_left = D_left, D_right = D_right))
  rsample <- rfun(n)
  
  # Perform Kolmogorov-Smirnov two-sample test 
  print(paste0("performing ks-test on ", fun_name,"..."))
  test <- ks.test(sample, rsample)
  p_value <- test$p.value
  error_msg <- paste0("test failed with ", fun_name)
  pass_msg <- paste0("test passed with ", fun_name)
  if(p_value <= 0.05)
    print(error_msg)
  else
    print(pass_msg)
}
  
##**************************************************************
## logistic distribution
dis_logistic <- function(x){
  return(exp(x)/(1+exp(x))^2)
}
test_fun(1000, dis_logistic, rlogis, "logistic distribution", D_left = -10, D_right = 10)


##***************************************************************
## standard normal distribution
dis_norm <- function(x) {
return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
}
test_fun(1000, dis_norm, rnorm, "normal distribution", D_left = -10, D_right = 10)


##***************************************************************
## uniform distribution
dis_unif <- function(x){
  if(x<=1&&x>=0) return(1)
  else
    return(0)
}
test_fun(1000, dis_unif, runif, "uniform distribution", D_left = 0, D_right = 1)


##***************************************************************
## laplace distribution
dis_laplace <- function(x, m = 0, s = 1) {
 return(exp(-abs(x-m)/s)/(2*s))
 }
test_fun(1000, dis_laplace, rmutil::rlaplace, "laplace distribution", D_left = -10, D_right = 10)


##***************************************************************
## gamma distribution
dis_gamma <- function(x, a = 1, b = 1){
  return((x^(a-1)*b^a*exp(-x*b))/gamma(a))
}
r_gamma <- function(n){
  return(rgamma(n, shape = 1))
}
test_fun(1000, dis_gamma, r_gamma, "gamma distribution", D_left = 0.01, D_right = Inf)


##***************************************************************
## non-concave case: a mixture of two normal distributions
dis_nc <- function(x){
  p <- 0.5*(1/sqrt(2*pi))*exp(-(x^2)/2) + 0.5*(1/sqrt(2*pi))*exp(-((x-3)^2)/2)
  return(p)
}
test_fun(1000, fun = dis_nc, rfun = NA , "mixture of normal distributions", D_left = -10, D_right = 10)








