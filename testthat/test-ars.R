#############
##  TESTS  ##
#############
  
context("tests")

##**************************************************************
## logistic distribution
test_that("logistic distribution",
          {
  test_fun(1000, dis_logistic, rlogis, "logistic distribution", D_left = -10, D_right = 10)
  expect_equal(1,1)
  })


##***************************************************************
## standard normal distribution
test_that("normal distribution",{
          test_fun(1000, dis_norm, rnorm, "normal distribution", D_left = -10, D_right = 10)
          expect_equal(1,1)
          })

##***************************************************************
## uniform distribution
test_that("uniform distribution",{
  test_fun(1000, dis_unif, runif, "uniform distribution", D_left = 0, D_right = 1)
  expect_equal(1,1)
  })

##***************************************************************
## laplace distribution
test_that("laplace distribution",{
  test_fun(1000, dis_laplace, rmutil::rlaplace, "laplace distribution", D_left = -10, D_right = 10)
  expect_equal(1,1)
  })

##***************************************************************
## gamma distribution
test_that("gamma distribution",{
  r_gamma <- function(n){
    return(rgamma(n, shape = 1))
  }
  test_fun(1000, dis_gamma, r_gamma, "gamma distribution", D_left = 0.01, D_right = Inf)
  expect_equal(1,1)
  })

##***************************************************************
## non-concave case: a mixture of two normal distributions
 # test_that("mixture of normal distribution",{
 # expect_error(test_fun(1000, fun = dis_nc, rfun = NA , "mixture of normal distributions", D_left = -10, D_right = 10))
 # })

