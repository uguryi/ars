#' Adaptive-Rejection Sampler Function
#'
#' This function allows you to perform a rejection sampling based on the 
#' theory developed by Gilks et al. (1992). This function employs the tangent
#' approach as opposed to the secant approach
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' ars()
#' 
#############
# LIBRARIES #
#############

#install.packages("assertthat")
# library(assertthat)

#install.packages("Deriv")
#library("Deriv")

#install.packages("distr")
# library(distr)
# 
# #install.packages("Ecfun")
# #library(Ecfun)
# 
# #install.packages("testthat")
# library(testthat)

################
# ARS FUNCTION #
################

ars <- function(n, 
                g = NA, 
                h = NA, 
                h_prime = NA, 
                D_left = -Inf, 
                D_right = Inf, 
                k = 10, 
                step = 3,
                center = 0) {
  
  ## Input:
  ## n: the number of desired samples 
  ## g: the density function of interest
  ## h: the log of g (users can provide either g or h)
  ## h_prime: h prime function (optional)
  ## D_left/D_right: the desired left/right end of domain (optional), default = -Inf/Inf
  ## k: the number of desired initial abscissae (optional), default = 10
  ## step: the desired width we use to choose left/right end if the domain is  
  ##       unbounded(optional), default = 3
  ## center: the desired center to search for left/right ends if the domain is  
  ##       unbounded(optional), default = 0
  
  #######################
  # INITIALIZATION STEP #
  #######################
  
  # Check the parameters
  assert_that(is.numeric(n), msg = "Please provide n as a number")
  assert_that(is.numeric(D_left), msg = "Please provide D_left as a number")
  assert_that(is.numeric(D_right), msg = "Please provide D_right as a number")
  assert_that(is.numeric(k), msg = "Please provide k as a number")
  assert_that(is.numeric(step), msg = "Please provide step as a number")
  assert_that(is.numeric(center), msg = "Please provide center as a number")
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
  
  
  # Initialize final sample, final sample counter, and update_needed flag
  final_sample <- numeric(n)
  count <- 1
  update_needed <- 0
  
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
      iter <- 1
      Check_boundary(D_right)
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
      Check_boundary(D_left)
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
    eps <- 1e-4
    flag <- all(h_prime_x[2:len] - h_prime_x[1:len-1] < eps)
    if(!flag) stop("Input is not a log-concave function", call. = FALSE)
  }
  check_concave(h_prime_x)
  
  # Calculate z
  calc_z <- function(k, D_left, D_right) {
    z <- numeric(k+1)
    z[1] <- ifelse(D_left == -Inf, T_k[1]-step, D_left)
    z[k+1] <- ifelse(D_right == Inf, T_k[k]+step, D_right)
    for (j in seq(2,k)) {
      # In case h is a linear function
      if(h_prime_x[j-1] - h_prime_x[j] < 1e-5) z[j] <- (T_k[j]+T_k[j-1])/2
      else
        z[j] <- (h_x[j] - h_x[j-1] - T_k[j] * h_prime_x[j] + T_k[j-1] * h_prime_x[j-1]) / 
          (h_prime_x[j-1] - h_prime_x[j])
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
  
  
  
  while (count <= n) {
    
    if (update_needed) {
      
      #################
      # UPDATING STEP #
      #################
      
      # Update T_k
      T_k <- c(T_k, x_star)
      T_k <- sort(T_k)
      
      # Update h_x and h_prime_x
      h_x_new <- h(T_k[l_interval+1])
      h_x <- append(h_x, h_x_new, after=l_interval)
      h_prime_x_new <- h_prime(T_k[l_interval+1])
      h_prime_x <- append(h_prime_x, h_prime_x_new, after=l_interval)
      
      # Update z
      calc_z_new <- function(l_interval) {
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
      z <- calc_z_new(l_interval)
      
      # Update denominator of s_j
      calc_denom_s_j_new <- function(l_interval) {
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
      denom <- calc_denom_s_j_new(l_interval)
      
      # Increment k
      k <- k + 1
      
    }
    
    #################
    # SAMPLING STEP #
    #################
    
    # Sample x_star from s_k
    piece_probs <- denom/sum(denom)
    piece_probs <- ifelse(piece_probs < 0, 0, piece_probs) # sometimes probs can be negative by mistake???
    piece_selected <- which(rmultinom(1, 1, piece_probs) != 0)
    x_star <- r(AbscontDistribution(d=calc_nom_s_j(piece_selected), 
                                    low1=z[piece_selected], 
                                    up1=z[piece_selected+1]))(1)
    
    # Determine where the sample x_star falls
    l_interval <- findInterval(x_star, T_k)
    u_interval <- findInterval(x_star, z)
    
    # Sample from Uniform(0,1)
    w <- runif(1)
    
    # Calculate l_k_x_star
    if (l_interval == 0 || l_interval == k) {
      l_k_x_star <- -Inf
    } else {
      l_k_x_star <- calc_l_j(l_interval)(x_star)
    }
    
    # Calculate u_k_x_star
    u_k_x_star <- calc_u_j(u_interval)(x_star)
    
    # Test whether to accept or reject the sample x_star
    if (w <= exp(l_k_x_star - u_k_x_star)) {
      final_sample[count] <- x_star
      count <- count + 1
      if (count <= n) {
        # Go back to sampling step
        next
      } else {
        return(final_sample)
      }
    } else {
      h_x_star <- h(x_star)
      h_prime_x_star <- h_prime(x_star)
      if (w <= exp(h_x_star - u_k_x_star)) {
        final_sample[count] <- x_star
        count <- count + 1
      }
      if (count <= n) {
        update_needed <- 1
        # Go back to sampling step
        next
      } else {
        return(final_sample)
      }
    }
  }
}
