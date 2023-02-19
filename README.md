# The PM Test 

This is the github repo for the PM test proposed by Fan et al. (2023) that tests overidentifying restrictions (or instrumental variable validity) with high-dimensional data and heteroskedastic errors. The baseline method is called the M test since it is based on estimation and inference for maximum norms of (high-dimensional) vectors. The recommended procedure for practitioners is called the power-enhanced M test (PM test), which enhances to power by a quadratic statistic. 

## Required packages 

You can install the required package with 

``` r
install.pacakges(c("glmnet,MASS,fastclime,matrixStats"))
```
## R Functions

"L1Penalty_functions.R"  contain necessary R functions for Lasso and CLIME estimators for bias correction. 

"example.R" provides an application example that is shown below.  

## Inputs and Outputs

```{r example}
# Inputs:
#    Y: n by 1 vector of the outcome variable.
#    D: n by 1 vector of the endogenous variable.
#    Z: n by pz matrix of instrumental variables.
#    X: n by px matrix of covariates. NULL by default.
# 
#    sig.level: Significance level.
#    intercept: Logical. FALSE for a model without intercept
#    tuning: Lasso tuning parameter selection method.
#    seed: random seed.
#    MM: Number of Bootstrapping replications.
#    A.mat: Weighting matrix A for the quadratic form of IV coefficients.
#    weight: A (p_x + p_z)-dimensional vector for Lasso estimation.
#    beta.output: If TRUE, export the point estimator and confidence interval of beta.
# 
# Outputs:
#   M_test: Logical. TRUE for the M test rejecting the validity of instruments.
#   PM_test: Logical. TRUE for the PM test rejecting the validity of instruments.
#   sig.level: Significance level. 0.05 by default.
#   M_p_val: P-value of the M test.
#   PM_p_val: P-value of the PM test.
#   beta.hat: (available only if bate.output == TRUE) Point estimator of beta.
#   beta.CI: (available only if bate.output == TRUE) (1-sig.level)*100% Confidence Interval of beta.
```





## Example

This is a basic example that shows you how to use the Q test. 

```{r example}
rm(list = ls())

library(glmnet)
library(fastclime)
library(matrixStats) 
library(MASS)

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
```


Test results with valid instruments 
```{r}
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

```


Test results with invalid instruments 
```{r}
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
```