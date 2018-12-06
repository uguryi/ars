###################
### CODE SKETCH ###
###################

# Author: Ugur Yildirim
# Date: Thursday, November 22, 2018

# Revised:
#  Tuesday, November 27, 2018
#  Wednesday, November 28, 2018
#  Thursday, November 29, 2018



#############
# LIBRARIES #
#############

#install.packages("Deriv")
library(Deriv)

#install.packages("Ecfun")
library(Ecfun)



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
init_T_k <- function(k, D_left, D_right, width) {
  T_k <- numeric(k)
  if (D_left == -Inf || D_right == Inf) {
    argmax <- optim(0, h, list(fnscale=-1))$par
  }
  if (D_left == -Inf) {
    T_k[1] <- argmax - width/2
  } else {
    T_k[1] <- D_left
  }
  if (D_right == Inf) {
    T_k[k] <- argmax + width/2
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
  nom <- rep(list(NA), k)
  denom <- numeric(k)
  for (j in seq(1:k)) {
    nom[[j]] <- function(x) {
      exp(h_x[j] + (x - T_k[j]) * h_prime_x[j])
    }
    denom[j] <- integrate(nom[[j]], lower = z[j], upper = z[j+1])
  }
  return(list(nom, denom))
}
nom <- calc_s_k(k)[[1]]
denom <- calc_s_k(k)[[2]]

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

while (count <= n) {
  
  # Sample x_star from s_k
  sample_from_s_k <- function() {
    piece_probs <- denom/sum(denom)
    piece_selected <- which(rmultinom(1, 1, piece_probs) != 0)
    lambda <- nom[[piece_selected]](0)
    x_star <- rtruncdist(1, 
                         rate = lambda, 
                         dist = 'exp', 
                         truncmin = z[piece_selected], 
                         truncmax = z[piece_selected + 1])
    return(x_star)
  } 
  x_star <- sample_from_s_k()
  
  # Determine where the sample x_star falls
  l_interval <- findInterval(x_star, T_k)
  u_interval <- findInterval(x_star, z)
  
  # Sample from Uniform(0,1)
  w <- runif(1)

  # Calculate l_k_x_star
  if (l_interval == 0 || l_interval == k) {
    l_k_x_star <- -Inf
  } else {
    l_k_x_star <- l_k[[l_interval]](x_star)
  }
  
  # Calculate u_k_x_star
  u_k_x_star <- u_k[[u_interval]](x_star)
  
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
          z <- append(z, (h_x[2] - h_x[1] - T_k[2] * h_prime_x[2] + T_k[1] * h_prime_x[1]) / (h_prime_x[1] - h_prime_x[2]), after=1)
          return(z)
        } else if (l_interval == k) {
          z <- append(z, (h_x[k+1] - h_x[k] - T_k[k+1] * h_prime_x[k+1] + T_k[k] * h_prime_x[k]) / (h_prime_x[k] - h_prime_x[k+1]), after=k)
          return(z)
        } else {
          z_new <- c()
          for (j in seq(l_interval+1:l_interval+2)) {
            z_new <- append(z_new, (h_x[j] - h_x[j-1] - T_k[j] * h_prime_x[j] + T_k[j-1] * h_prime_x[j-1]) / (h_prime_x[j-1] - h_prime_x[j]))
          }
          z[l_interval+1] <- z_new[1]
          z <- append(z, z_new[2], after=l_interval+1)
          return(z)
        }
      }
      z <- calc_z_new(l_interval)
      
      # Update u_k
      calc_u_k_new <- function(l_interval) {
        if (l_interval == 0) {
          u_k_new <- list()
          for (j in seq(1:2)) {
            u_k_new <- append(u_k_new, function(x) h_x[j] + (x - T_k[j]) * h_prime_x[j])
          }
          u_k[1] <- u_k_new[[2]]
          u_k <- append(u_k, u_k_new[[1]], after=0)
          return(u_k)
        } else if (l_interval == k) {
          u_k_new <- list()
          for (j in seq(k:k+1)) {
            u_k_new <- append(u_k_new, function(x) h_x[j] + (x - T_k[j]) * h_prime_x[j])
          }
          u_k[k] <- u_k_new[[1]]
          u_k <- append(u_k, u_k_new[[2]])
          return(u_k)
        } else {
          u_k_new <- list()
          for (j in seq(l_interval:l_interval+2)) {
            u_k_new <- append(u_k_new, function(x) h_x[j] + (x - T_k[j]) * h_prime_x[j])
          }
          u_k[l_interval] <- u_k_new[[1]]
          u_k[l_interval+1] <- u_k_new[[3]]
          u_k <- append(u_k, u_k_new[[2]], after=l_interval)
          return(u_k)
        }
      }
      u_k <- calc_u_k_new(l_interval)
      
      # Update s_k
      calc_s_k_new <- function(l_interval) {
        if (l_interval == 0) {
          nom_new <- list()
          denom_new <- numeric(2)
          for (j in seq(1:2)) {
            nom_new <- append(nom_new, function(x) exp(h_x[j] + (x - T_k[j]) * h_prime_x[j]))
          }
          denom_new[1] <- integrate(nom[[1]], lower = z[1], upper = z[2])
          denom_new[2] <- integrate(nom[[2]], lower = z[2], upper = z[3])
          nom[1] <- nom_new[[2]]
          nom <- append(nom, nom_new[[1]], after=0)
          denom[1] <- denom_new[2]
          denom <- append(denom, denom_new[1], after=0)
          return(list(nom, denom))
        } else if (l_interval == k) {
          nom_new <- list()
          denom_new <- numeric(2)
          for (j in seq(k:k+1)) {
            nom_new <- append(nom_new, function(x) exp(h_x[j] + (x - T_k[j]) * h_prime_x[j]))
          }
          denom_new[1] <- integrate(nom[[1]], lower = z[k], upper = z[k+1])
          denom_new[2] <- integrate(nom[[2]], lower = z[k+1], upper = z[k+2])
          nom[k] <- nom_new[[1]]
          nom <- append(nom, nom_new[[2]])
          denom[k] <- denom_new[1]
          denom <- append(denom, denom_new[2])
          return(list(nom, denom))
        } else {
          nom_new <- list()
          denom_new <- numeric(3)
          for (j in seq(l_interval:l_interval+2)) {
            nom_new <- append(nom_new, function(x) exp(h_x[j] + (x - T_k[j]) * h_prime_x[j]))
          }
          denom_new[1] <- integrate(nom[[1]], lower = z[l_interval], upper = z[l_interval+1])
          denom_new[2] <- integrate(nom[[2]], lower = z[l_interval+1], upper = z[l_interval+2])
          denom_new[3] <- integrate(nom[[3]], lower = z[l_interval+2], upper = z[l_interval+3])
          nom[l_interval] <- nom_new[[1]]
          nom[l_interval+1] <- nom_new[[3]]
          nom <- append(nom, nom_new[[2]], after=l_interval)
          denom[l_interval] <- denom_new[1]
          denom[l_interval+1] <- denom_new[3]
          denom <- append(denom, denom_new[2], after=l_interval)
          return(list(nom, denom))
        }
      }
      nom <- calc_s_k_new(l_interval)[[1]]
      denom <- calc_s_k_new(l_interval)[[2]]
      
      # Update l_k
      calc_l_k_new <- function(l_interval) {
        if (l_interval == 0) {
          l_k <- append(l_k, ((T_k[2] - x) * h_x(T_k[1]) + (x - T_k[1]) * h_x(T_k[2])) / (T_k[2] - T_k[1]), after=0)
          return(l_k)
        } else if (l_interval == k) {
          l_k <- append(l_k, ((T_k[k+1] - x) * h_x(T_k[k]) + (x - T_k[k]) * h_x(T_k[k+1])) / (T_k[k+1] - T_k[k]))
          return(l_k)
        } else {
          l_k_new <- list()
          for (j in seq(l_interval:l_interval+1)) {
            l_k_new <- append(l_k_new, function(x) ((T_k[j+1] - x) * h_x(T_k[j]) + (x - T_k[j]) * h_x(T_k[j+1])) / (T_k[j+1] - T_k[j]))
          }
          l_k[l_interval] <- l_k_new[1]
          l_k <- append(l_k, l_k_new[2], after=l_interval)
          return(l_k)
        }
      }
      l_k <- calc_l_k_new(k)
      
      # Increment k
      k <- k + 1
      
      # Go back to sampling step
      next
    } else {
      return(final_sample)
    }
  }
}
