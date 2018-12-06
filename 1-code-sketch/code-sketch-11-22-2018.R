###################
### CODE SKETCH ###
###################

# Author: Ugur Yildirim
# Date: Thursday, November 22, 2018





#############
# LIBRARIES #
#############

#install.packages("Deriv")
library(Deriv)





#######################
# INITIALIZATION STEP #
#######################

# Initialize final sample and final sample counter

n <- ... # number of samples we want
final_sample <- numeric(n)
count <- 1



# Define h and h'

h <- function(x) {...} # whatever the function h is
h_prime <- Deriv(h)



# Initialize T_k

k <- ... # number of x_i we want
D_left <- ... # left end of domain
D_right <- ... # right end of domain

init_T_k <- function(k, D_left, D_right) {
  
  T_k <- numeric(k)
  
  if (D_left == -Inf) {
    T_k[1] <- ... # some value s.t. h_prime(T_k[1]) > 0 
    # WE NEED TO FIGURE OUT HOW TO DO THIS EXACTLY
  } else {
    T_k[1] <- D_left
  }
  
  if (D_right == Inf) {
    T_k[k] <- ... # some value s.t. h_prime(T_k[k]) < 0 
    # WE NEED TO FIGURE OUT HOW TO DO THIS EXACTLY
  } else {
    T_k[k] <- D_right
  }
  
  step <- (T_k[k] - T_k[1])/(k-1)
  
  for (i in seq(2:(k-1))) {
    T_k[i] <- T_k[i-1] + step
  }
  
  return(T_k)
  
}

T_k <- init_T_k(k, D_left, D_right)



# Evaluate h and h' at T_k

h_x <- h(T_k)
h_prime_x <- h_prime(T_k)



# Calculate z

calc_z <- function(k, D_left, D_right) {
  
  z <- numeric(k+1)
  
  z[1] <- D_left
  z[k+1] <- D_right
  
  for (j in seq(2:k)) {
    
    z[j] <- (h_x[j] - h_x[j-1] - T_k[j] * h_prime_x[j] + T_k[j-1] * h_prime_x[j-1]) / (h_prime_x[j-1] - h_prime_x[j])
    
  }
  
  return(z)
  
}

z <- calc_z(k, D_left, D_right)



# Calculate u_k

calc_u_k <- function(k) {
  
  u_k <- rep(list(NA), k)
  
  for (j in seq(1:k)) {
    
    u_k[[j]] <- function(x) {
      
      h_x[j] + (x - T_k[j]) * h_prime_x[j]
      
    }
    
  }
  
  return(u_k)
  
}

u_k <- calc_u_k(k)



# Calculate s_k

calc_s_k <- function(k) {
  
  s_k <- rep(list(NA), k)
  
  for (j in seq(1:k)) {
    
    nom <- function(x) {
      
      exp(h_x[j] + (x - T_k[j]) * h_prime_x[j])
      
    }
    
    denom <- integrate(nom, lower = z[i], upper = z[i+1])
    
    s_k[[j]] <- (exp(h_x[j] + (x - T_k[j]) * h_prime_x[j])) / denom
    # NOT SURE IF IT'S VALID TO DO THIS ON A PIECE BY PIECE BASIS LIKE THIS
    
  }
  
  return(s_k)
  
}

s_k <- calc_s_k(k)



# Calculate l_k

calc_l_k <- function(k) {
  
  l_k <- rep(list(NA), k-1)
  
  for (j in seq(1:(k-1))) {
    
    l_k[[j]] <- function(x) {
      
      ((T_k[j+1] - x) * h_x(T_k[j]) + (x - T_k[j]) * h_x(T_k[j+1])) / (T_k[j+1] - T_k[j])
      
    } 
    
  }
  
  return(l_k)
  
}

l_k <- calc_l_k(k)





#################
# SAMPLING STEP #
#################

sample_from_s_k <- function() {...} 
# I DON'T KNOW HOW TO SAMPLE FROM A PIECEWISE EXPONENTIAL FUNCTION...

x_star <- sample_from_s_k()

l_interval <- findInterval(x_star, T_k)
u_interval <- findInterval(x_star, z)

w <- runif(1)

test_x_star <- function(l_interval, u_interval) {
  
  if (l_interval == 0 || l_interval == k) {
    
    l_k_x_star <- -Inf
    
  } else {
    
    l_k_x_star <- l_k[[l_interval]](x_star)
    
  }
  
  u_k_x_star <- u_k[[u_interval]](x_star)
  
  if (w <= exp(l_k_x_star - u_k_x_star)) {
    
    final_sample[count] <- x_star
    count <- count + 1
    
    if (count <= n) {
      
      # Go back to Sampling Step
      
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
      
      # Go to Updating Step
      
    } else {
      
      return(final_sample)
      
    }
    
  }
  
}





#################
# UPDATING STEP #
#################

T_k <- c(T_k, x_star)

T_k <- sort(T_k)

h_x <- h(T_k)
h_prime_x <- h_prime(T_k)
# WE ACTUALLY SHOULDN'T HAVE TO EVALUATE THESE OVER AND OVER AGAIN...

z <- calc_z(k+1, D_left, D_right)

u_k <- calc_u_k(k+1)
s_k <- calc_s_k(k+1)
l_k <- calc_l_k(k+1)

k <- k + 1

# Go back to Sampling Step
