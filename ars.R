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
  
  
  #source("auxiliary.R")
  
  #######################
  # INITIALIZATION STEP #
  #######################
  
  # Check the parameters
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
  denom <- sapply(1:k, calc_denom_s_j)
 
  
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
