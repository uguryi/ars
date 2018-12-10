#' Adaptive-Rejection Sampling
#'
#' The objective of this package is to create a function that
#' performs adaptive-rejection sampling based on the tangent approach described 
#' by Gilks et al. (1992).
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
#' @return Results from the adaptive rejection sample from the information provided by the user
#' 
#' @examples 
#' ars(n=1000, g=dnorm , D_left = -Inf, D_right = Inf, k = 20)
#' 
#' @export
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
  
  
  #source("auxiliary.R")
  
  #######################
  # INITIALIZATION STEP #
  #######################
  
 # Check the parameters
  check_param(n, g, h, h_prime, D_left, D_right, k, step, center)
  
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
  if(identical(h_prime, NA)){
    h_prime <- function(x){
      if(length(x) > 1)  return(sapply(x,FUN = prime, fun = h, a = D_left, b = D_right))
      return(prime(x,h,D_left,D_right))
    }
  }
  
  # Initialize final sample, final sample counter, and update_needed flag
  final_sample <- numeric(n)
  count <- 1
  update_needed <- 0
  
  # Initialize T_k
  T_k <- init_T_k(k, D_left, D_right, h_prime, step, center)
  
  # Evaluate h and h' at T_k
  h_x <- sapply(T_k,h)
  h_prime_x <- h_prime(T_k)
  
  # Check concavity of h
  check_concave(h_prime_x)
  
  # Calculate z
  z <- calc_z(k, h_x, h_prime_x, T_k, D_left, D_right, step)
  
  # Calculate denominator of s_j
  denom <- sapply(1:k, calc_denom_s_j, h_x, T_k, h_prime_x, z)
  
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
      
      # Check concavity for updated h_prime
      check_concave(h_prime_x)
      
      # Update z
      z <- calc_z_new(z, k, l_interval, h_x, h_prime_x, T_k)
      
      # Update denominator of s_j
      denom <- calc_denom_s_j_new(z, k, denom, l_interval, h_x, h_prime_x, T_k)
      
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
    x_star <- distr::r(distr::AbscontDistribution(d=calc_nom_s_j(piece_selected, h_x, T_k, h_prime_x), 
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
      l_k_x_star <- calc_l_j(l_interval, T_k, h_x)(x_star)
    }
    
    # Calculate u_k_x_star
    u_k_x_star <- calc_u_j(u_interval, h_x, T_k, h_prime_x)(x_star)
    
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
