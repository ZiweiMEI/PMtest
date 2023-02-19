rm(list = ls())

library(glmnet)
library(fastclime)
library(matrixStats) 

source("PMtest.R")

n = 300
px = 250
pz = 100 

p = px + pz


phi <- c(0.5^(0:9),rep(0,px-10))
psi <- c(0.6^(0:9),rep(0,px-10))

s1 = 7 # num of relevant IV
gamma <- matrix(c(rep(1,s1),rep(0,pz-s1)),ncol = 1)

beta = 1

set.seed(2023)

## covariates and instruments 
Sigma = diag(p)
W <- MASS::mvrnorm(n=n, rep(0, p), Sigma)
if (px == 0){
  X <- NULL
  Z <- W
}else{
  X <- W[,1:px]
  Z <- W[,(px+1):p]
}

## error terms 
err <- MASS::mvrnorm(n=n, rep(0, 2), matrix( c(1, .5, .5,1),2))
e <- err[,1]
eps2 <- err[,2]

## Test result with valid instruments 

D <-  X%*%psi  + Z %*% gamma + eps2
Y <-   D*beta + X%*%phi + e

PMtest(Y,D,Z,X) 

# $M_test
# [1] FALSE
# 
# $PM_est
# [1] FALSE
# 
# $sig.level
# [1] 0.05
# 
# $p_val_M
# [1] 0.6233
# 
# $p_val_PM
# [1] 0.6233
# 
# $beta.hat
# [1] 0.9933557
# 
# $CI
# [1] 0.950848 1.035863

## Test result with invalid instruments 

s2 = 3 # number of invalid IVs
pi <-  c(rep(0.3,s2),rep(0,pz-s2))
Y <- D*beta + X%*%phi + Z %*% pi + e

PMtest(Y,D,Z,X) 

# $M_test
# [1] TRUE
# 
# $PM_est
# [1] TRUE
# 
# $sig.level
# [1] 0.05
# 
# $p_val_M
# [1] 0.013
# 
# $p_val_PM
# [1] 0
# 
# $beta.hat
# [1] 1.129874
# 
# $CI
# [1] 1.091233 1.168515