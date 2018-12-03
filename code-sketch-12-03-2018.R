#######################
# INITIALIZATION STEP #
#######################

# Shubei Wang
# Dec 03, 2018

## Input:
## n: the number of desired samples 
## g: the density function of interest
## h: the log of g (users can provide either g or h)
## h_prime: h prime function (optional)
## D_left/D_right: the desired left/right end of domain (optional), default = -Inf/Inf
## k: the number of desired initial abscissae (optional), default = 100
## step: the desired width we use to choose left/right end if the domain is  
##       unbounded(optional), default = 1

#############
# LIBRARIES #
#############

#install.packages("Ecfun")
library(Ecfun)

#install.packages("assertthat")
library(assertthat)

ars <- function(n,
                g = NA,
                h = NA,
                h_prime = NA,
                D_left = -Inf,
                D_right = Inf,
                k = 100,
                step = 1
){
  # Check the parameters
  assert_that(is.numeric(n), msg = "Please provide n as a number")
  assert_that(is.numeric(k), msg = "Please provide k as a number")
  assert_that(is.numeric(step), msg = "Please provide step as a number")
  assert_that(is.numeric(D_left), msg = "Please provide D_left as a number")
  assert_that(is.numeric(D_right), msg = "Please provide D_right as a number")
  if(D_left==D_right)
    stop("Please provide different D_left, D_right", call. = FALSE)
  if(identical(g, NA)&&identical(h,NA))
    stop("Please provide either g or h", call. = FALSE)
  if(!identical(g, NA)){
    if(class(g) != "function"){
      stop("Please provide g as a function", call. = FALSE)
    }
  }
  if(!identical(h, NA)){
    if(class(h) != "function"){
      stop("Please provide h as a function", call. = FALSE)
    }
  }
  if(!identical(h_prime, NA)){
    if(class(h_prime) != "function"){
      stop("Please provide h_prime as a function", call. = FALSE)
    }
  }
  
  
  # If h is not provided, calculate the log of g 
  if(identical(h, NA)){
    h <- function(x){
      return(log(g(x)))
    }
  }
  
  
  # If h_prime is not provided, calculate it numerically
  
  # x is the point where derivative is evaluated
  # fun is the function of interest
  # a, b are the left/right ends of domain
  prime <- function(x, fun, a, b){
    if(x == a) return((fun(x+1e-10)-fun(x))/1e-10)
    if(x == b) return((fun(x)-fun(x-1e-10))/1e-10)
    if(x>a && x<b) return((fun(x+1e-10)-fun(x-1e-10))/2e-10)
  }
  if(identical(h_prime,NA)){
    h_prime <- function(x){
      if(length(x)>1)  return(sapply(x,FUN = prime, fun = h, a = D_left, b = D_right))
      return(prime(x,h,D_left,D_right))
    }
  }
  
  
  # Initialize final sample and final sample counter
  final_sample <- numeric(n)
  count <- 1
  
  # Check if h is defined on boundaries
  Check_boundary <- function(x){
    if(!is.finite(h_prime(x))) {
      stop('Boundary is not defined on h', call. = FALSE)
    }
  }
  
  # Initialize T_k
  init_T_k <- function(k, D_left, D_right) {
    T_k <- numeric(k)
    
    # Check if h is defined on boundaries
    if(D_left != -Inf && D_right != Inf){
      Check_boundary(D_left)
      Check_boundary(D_right)
    }
    
    # If the left end is unbounded
    if(D_left == -Inf && D_right != Inf){
      Check_boundary(D_right)
      D_left <- D_right-step
      signal <- h_prime(D_left)
      while(signal <= 0){
        D_left <- D_left-step
        signal <- h_prime(D_left)
      }
    }
    
    # If the right end is unbounded
    if(D_left != -Inf && D_right == Inf){
      Check_boundary(D_left)
      D_right <- D_left+step
      signal <- h_prime(D_right)
      while(signal >= 0){
        D_right <- D_right+step
        signal <- h_prime(D_right)
      }
    }
    
    # If both are unbounded
    if(D_left == -Inf && D_right == Inf){
      start <- 0
      D_left <- start-step
      D_right <- start+step
      signal1 <- h_prime(D_left)
      signal2 <- h_prime(D_right)
      while((signal1 <= 0 || signal2 >= 0)&&is.finite(signal1)&&is.finite(signal2)){
        D_left <- D_legt-step
        D_right <- D_right+step
        signal1 <- h_prime(D_left)
        signal2 <- h_prime(D_right)
      }
    }
     
    T_k[1] <- D_left
    T_k[k] <- D_right
    len <- (T_k[k] - T_k[1])/(k-1)
    for (i in seq(2,(k-1))) {
      T_k[i] <- T_k[i-1] + len
    } 
    return(T_k)
  }
  T_k <- init_T_k(k, D_left, D_right)
  
  
  
  # Evaluate h and h' at T_k
  h_x <- sapply(T_k,h)
  h_prime_x <- h_prime(T_k)
  
  
  # Check concavity of h
  check_concave <- function(h_prime_x){
    len <- length(h_prime_x)
    flag <- all(h_prime_x[2:len] <= h_prime_x[1:len-1])
    if(!flag) stop("Input is not a log-concave function", call. = FALSE)
  }
  
  
  # Calculate z
  calc_z <- function(k, D_left, D_right) {
    z <- numeric(k+1)
    z[1] <- D_left
    z[k+1] <- D_right
    for (j in seq(2,k)) {
      # In case h is a linear function
      if(h_prime_x[j-1] - h_prime_x[j] < 1e-5) z[j] <- (T_k[j]+T_k[j-1])/2
      else
        z[j] <- (h_x[j] - h_x[j-1] - T_k[j] * h_prime_x[j] + T_k[j-1] * h_prime_x[j-1]) / (h_prime_x[j-1] - h_prime_x[j])
    }
    return(z)
  }
  z <- calc_z(k, D_left, D_right)
  
  # Define function to calculate u_j
  calc_u_j <- function(j) {
    u_j <- function(x) {
      h_x[j] + (x - T_k[j]) * h_prime_x[j]
    }
    return(u_j)
  }
  
  # Define function to calculate denominator of s_j
  calc_denom_s_j <- function(j) {
    s_j <- function(x) {
      exp(h_x[j] + (x - T_k[j]) * h_prime_x[j])
    }
    denom_s_j <- integrate(s_j, lower = z[j], upper = z[j+1])$value
    return(denom_s_j)
  }
  denom <- numeric(k)
  for (j in seq(1,k)) {
    denom[j] <- calc_denom_s_j(j)
  }
  
  # Define function to calculate nominator of s_j
  calc_nom_s_j <- function(j) {
    nom_s_j <- function(x) {
      exp(h_x[j] + (x - T_k[j]) * h_prime_x[j])
    }
    return(nom_s_j)
  }
  
  # Define function to calculate l_j
  calc_l_j <- function(j) {
    l_j <- function(x) {
      ((T_k[j+1] - x) * h_x[j] + (x - T_k[j]) * h_x[j+1]) / (T_k[j+1] - T_k[j])
    } 
    return(l_j)
  }
}

##****************************************************************
## I tested the code above using standard normal:
## dnorm <- function(x) {
## return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
## }
## ars(n=100, g=dnorm , D_left = -1, D_right = 1)
## ars(n=100, g=dnorm )
## 
## and uniform case:
## unif <- function(x){
## if(x<=1&&x>=0) return(1)
## else
##  return(0)
## }
## ars(n=100,g=unif,D_left = 0,D_right = 1)
##
## and laplace distribution:
## dlaplace <- function(x, m = 0, s = 1) {
## return(exp(-abs(x-m)/s)/(2*s))
## }
## ars(n=100,g=dlaplace)
##*******************************************************************



