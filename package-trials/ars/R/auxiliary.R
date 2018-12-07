#' Adaptive-Rejection Sampling Parameters Test
#'
#' The objective of this function is to test the parameter inputs of the user. It checks
#' for validity of the inputs and returns comments if an error is encountered with the inputs.
#'
#' @param n The number of desired samples
#' @param g The density function of interest
#' @param h The log of g (users can provide either g or h)
#' @param h_prime H prime function (optional)
#' @param D_left The desired left end of domain (optional), default = -Inf
#' @param D_right The desired right end of domain (optional), default = Inf
#' @param k The number of desired initial abscissae (optional), default = 10
#' @param step The desired width (optional), default = 3
#' @param center the desired center (optional), default = 0
#' 
#' @return If an error is encounter, provides the user some feedback regarding the source of error.
#'
#' @examples 
#'
#'


######################
# AUXILIARY FUNCTION #
######################

# Define the function to check parameters
check_param <- function(n, 
                        g,
                        h, 
                        h_prime, 
                        D_left, 
                        D_right, 
                        k, 
                        step,
                        center){
  
assertthat::assert_that(is.numeric(n), msg = "Please provide n as a number")
assertthat::assert_that(is.numeric(D_left), msg = "Please provide D_left as a number")
assertthat::assert_that(is.numeric(D_right), msg = "Please provide D_right as a number")
assertthat::assert_that(is.numeric(k), msg = "Please provide k as a number")
assertthat::assert_that(is.numeric(step), msg = "Please provide step as a number")
assertthat::assert_that(is.numeric(center), msg = "Please provide center as a number")
if(D_left == D_right)
  stop("Please provide different D_left, D_right", call. = FALSE)
if(identical(g, NA) && identical(h,NA))
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
}}

# Define function to initialize T_k
init_T_k <- function(k, D_left, D_right, h_prime, step, center) {

  # Check if h is defined on boundaries
  if(D_left != -Inf && D_right != Inf){
    
    if(!is.finite(h_prime(D_left))) {
      stop('Left boundary is not defined on h', call. = FALSE)
    }
    if(!is.finite(h_prime(D_right))) {
      stop('Right boundary is not defined on h', call. = FALSE)
    }
  }
  
  # If the left end is unbounded
  if(D_left == -Inf && D_right != Inf){
    iter <- 1
    if(!is.finite(h_prime(D_right))) {
      stop('Right boundary is not defined on h', call. = FALSE)
    }
    D_left <- D_right-step
    signal <- h_prime(D_left)
    while(signal <= 0 && iter < 50){
      D_left <- D_left-step
      signal <- h_prime(D_left)
      iter <- iter+1
    }
    if(iter == 50)
      stop("Please provide valid boundaries", call. = FALSE)
  }
  
  # If the right end is unbounded
  if(D_left != -Inf && D_right == Inf){
    iter <- 1
    if(!is.finite(h_prime(D_left))) {
      stop('Left boundary is not defined on h', call. = FALSE)
    }
    D_right <- D_left+step
    signal <- h_prime(D_right)
    while(signal >= 0 && iter < 50){
      D_right <- D_right+step
      signal <- h_prime(D_right)
      iter <- iter+1
    }
    if(iter == 50)
      stop("Please provide valid boundaries", call. = FALSE)
  }
  
  # If both are unbounded
  if(D_left == -Inf && D_right == Inf){
    iter <- 1
    D_left <- center-step
    D_right <- center+step
    signal1 <- h_prime(D_left)
    signal2 <- h_prime(D_right)
    while((signal1 <= 0 || signal2 >= 0)&& iter < 50){
      D_left <- D_legt-step
      D_right <- D_right+step
      signal1 <- h_prime(D_left)
      signal2 <- h_prime(D_right)
      iter <- iter+1
    }
    if(iter == 50)
      stop("Please provide valid center", call. = FALSE)
  }
  
  T_k <- seq(D_left, D_right, length.out =  k)
  return(T_k)
}

# Define function to check concavity of h
check_concave <- function(h_prime_x){
  len <- length(h_prime_x)
  eps <- 1e-4
  flag <- all(h_prime_x[2:len] - h_prime_x[1:len-1] < eps)
  if(!flag) stop("Input is not a log-concave function", call. = FALSE)
}

# Define function to calculate z
calc_z <- function(k, h_x, h_prime_x, T_k, D_left, D_right, step){
  
  zl <- ifelse(D_left == -Inf, T_k[1]-step, D_left)
  zr <- ifelse(D_right == Inf, T_k[k]+step, D_right)
  
  h1 <- h_x[1:k-1]
  h2 <- h_x[2:k]
  hp1 <- h_prime_x[1:k-1]
  hp2 <- h_prime_x[2:k]
  t1 <- T_k[1:k-1]
  t2 <- T_k[2:k]
  
  check <- 1/(hp1 - hp2)
  check[!is.finite(check)] <- 0
  z <- ((hp1 - hp2) < 1e-5)*0.5*(t1+t2) + ((hp1 - hp2) >= 1e-5)*(h2 - h1 - t2*hp2 + t1*hp1)*check
  z<- c(zl, z, zr)
  return(z)
}

# Define function to update z
calc_z_new <- function(z, k, l_interval, h_x, h_prime_x, T_k) {
  if (l_interval == 0) {
    
    # In case h is a linear function
    if(h_prime_x[1] - h_prime_x[2] < 1e-5) 
      z_new <- (T_k[1]+T_k[2])/2
    else
      z_new <- (h_x[2] - h_x[1] - T_k[2] * h_prime_x[2] + T_k[1] * h_prime_x[1])/(h_prime_x[1] - h_prime_x[2])
    z <- append(z, z_new, after=1)
    return(z)
  } else if (l_interval == k) {
    
    # In case h is a linear function
    if(h_prime_x[k] - h_prime_x[k+1] < 1e-5) 
      z_new <- (T_k[k]+T_k[k+1])/2
    else
      z_new <- (h_x[k+1] - h_x[k] - T_k[k+1] * h_prime_x[k+1] + T_k[k] * h_prime_x[k])/(h_prime_x[k] - h_prime_x[k+1])
    z <- append(z, z_new, after=k)
    return(z)
  } else {
    
    z_new <- c()
    for (j in seq(l_interval+1,l_interval+2)) {
      
      # In case h is a linear function
      if(h_prime_x[j-1] - h_prime_x[j] < 1e-5) 
        new <- (T_k[j]+T_k[j-1])/2
      else
        new <- (h_x[j] - h_x[j-1] - T_k[j] * h_prime_x[j] + T_k[j-1] * h_prime_x[j-1])/(h_prime_x[j-1] - h_prime_x[j])
      
      z_new <- append(z_new, new)
    }
    z[l_interval+1] <- z_new[1]
    z <- append(z, z_new[2], after=l_interval+1)
    return(z)
  }
}

# Update denominator of s_j
calc_denom_s_j_new <- function(z, k, denom, l_interval, h_x, h_prime_x, T_k) {
  if (l_interval == 0) {
    denom_new <- numeric(2)
    for (j in seq(1,2)) {
      denom_new[j] <- integrate(function(x) exp(h_x[j] + (x - T_k[j]) * h_prime_x[j]), 
                                lower = z[j], 
                                upper = z[j+1])$value
    }
    denom[1] <- denom_new[2]
    denom <- append(denom, denom_new[1], after=0)
    return(denom)
  } else if (l_interval == k) {
    denom_new <- numeric(2)
    for (j in seq(k,k+1)) {
      denom_new[j-k+1] <- integrate(function(x) exp(h_x[j] + (x - T_k[j]) * h_prime_x[j]), 
                                    lower = z[j], 
                                    upper = z[j+1])$value
    }
    denom[k] <- denom_new[1]
    denom <- append(denom, denom_new[2])
    return(denom)
  } else {
    denom_new <- numeric(3)
    denom_new[1] <- integrate(function(x) exp(h_x[l_interval] + (x - T_k[l_interval]) * h_prime_x[l_interval]), 
                              lower = z[l_interval], 
                              upper = z[l_interval+1])$value
    denom_new[2] <- integrate(function(x) exp(h_x[l_interval+1] + (x - T_k[l_interval+1]) * h_prime_x[l_interval+1]), 
                              lower = z[l_interval+1], 
                              upper = z[l_interval+2])$value
    denom_new[3] <- integrate(function(x) exp(h_x[l_interval+2] + (x - T_k[l_interval+2]) * h_prime_x[l_interval+2]), 
                              lower = z[l_interval+2], 
                              upper = z[l_interval+3])$value
    denom[l_interval] <- denom_new[1]
    denom[l_interval+1] <- denom_new[3]
    denom <- append(denom, denom_new[2], after=l_interval)
    return(denom)
  }
}
















