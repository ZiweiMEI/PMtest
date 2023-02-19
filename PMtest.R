# Packages required: 
#   glmnet, MASS, fastclime, matrixStats, rockchalk
# 
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
#    rockchalk: TRUE for using the rockchalk package for replication of the results from the mvrnorm function. 
#    
# Outputs: 
#   M_test: Logical. TRUE for the M test rejecting the validity of instruments. 
#   PM_test: Logical. TRUE for the PM test rejecting the validity of instruments.
#   sig.level: Significance level. 0.05 by default.
#   M_p_val: P-value of the M test. 
#   PM_p_val: P-value of the PM test. 
#   beta.hat: (available only if bate.output == TRUE) Point estimator of beta.
#   beta.CI: (available only if bate.output == TRUE) (1-sig.level)*100% Confidence Interval of beta.


PMtest <- function(Y, D, Z, X = NULL, 
                  sig.level = 0.05,
                  intercept = F,
                  tuning = "CV", 
                  seed = NULL,
                  MM = 10000, 
                  A.mat = NULL,
                  weight = NULL,
                  beta.output = TRUE){
  
  
  
  
  if (intercept){
    if (is.null(X)){
      X = matrix(1,n,1)
    }else{
    X = cbind(1,X)
    }
    intercept = F 
  }
  
  if (is.null(X)){
    px = 0
  }else{
    px = dim(X)[2]
  }
  
  W = cbind(X,Z)
  
  if (is.null(weight)){
    weight = rep(1,ncol(W))
  }
  
  n = dim(Z)[1]
  pz = dim(Z)[2]
  p = px + pz 
  
  
  Sigma.hat = t(W)%*%W/n
  
  Sigma.hat <- 1/2 * (Sigma.hat + t(Sigma.hat))
  
   
  ############# Auxilliary parameter estimation 
  ## estimate \gamma 
  B.hat <-  Lasso(W,D,standardize = TRUE,intercept = intercept,tuning = tuning,weight = weight)
  e2_h = D - W %*% B.hat 
  gamma.hat <- B.hat[(px+1):p,] 
  ## \hat{u}_2 
  SigmaZ.hat = diag(diag(t(Z)%*%Z/n))
  
  if (is.null(A.mat)){
    A.mat = SigmaZ.hat
  }
  out1 = fastclime(W)
  out2 = fastclime.selector(out1$lambdamtx, out1$icovlist, min(out1$lambdamtx[out1$lambdamtx>0]))
  Omega.hat = out2$icov
  u2 = Omega.hat %*% matrix(c(rep(0,px),A.mat %*% gamma.hat),ncol = 1)
  
  ## estimate \|\gamma\|_2^2
  gamma.norm.sq.hat <- as.numeric(t(gamma.hat) %*% A.mat %*% gamma.hat) + 2/n * sum((W %*% u2) * (D-W%*%B.hat))
  
  
  ## estimate \Gamma 
  A.hat <- Lasso(W,Y,standardize = TRUE,intercept = intercept,tuning = tuning,weight = weight)
  
  e1_h = Y - W %*% A.hat 
  Gamma.hat <- A.hat[(px+1):p,]
 
  
  # \hat{u}_1  
  u1 = Omega.hat %*% matrix(c(rep(0,px),A.mat %*% Gamma.hat),ncol = 1)
  
  ## estimate <\Gamma,\gamma>
  GammaTgamma.hat <- as.numeric(t(gamma.hat) %*% A.mat %*% Gamma.hat) + 1/n * sum((W %*% u1) * (D-W%*%B.hat)) +
    1/n * sum((W %*% u2) * (Y-W%*%A.hat))
  
  ## estimate \beta 
  if (gamma.norm.sq.hat <= 0){
    beta.hat = 0
  }else{beta.hat <- GammaTgamma.hat / gamma.norm.sq.hat}
  
  e1.hat <- Y - W %*% A.hat
  e2.hat <- D - W %*% B.hat
  e.hat <- e1.hat - beta.hat * e2.hat
  
  if (gamma.norm.sq.hat > 0){
    sd.beta =  sqrt(1/n * sum( as.numeric(W %*% u2)^2*(e.hat)^2) ) / gamma.norm.sq.hat
  }else{
    sd.beta = 0
  } 
 
  CI.up = beta.hat + sd.beta * qnorm(1 - sig.level/2 ) / sqrt(n)
  CI.low = beta.hat - sd.beta * qnorm(1 - sig.level/2 )  / sqrt(n)
  
  
  ############## Bias-corrected pi 
  
  
  zeta.hat <- Y - beta.hat * D
   
  C.hat <- Lasso(W,zeta.hat,standardize = TRUE,intercept = intercept,tuning = tuning, weight = weight) 
  pi.hat <- C.hat[(px+1):p,]      
  
  
  
  C.hat.db = C.hat +  Omega.hat %*% ((1/n) * t(W) %*% (zeta.hat  - W %*% C.hat )) 
  pi.hat.db = C.hat.db[(px+1):p,]
  
  
  e_h <- as.numeric(zeta.hat - W %*% C.hat)
  Sigmae = crossprod( (W*e_h) ) / n 
  if (gamma.norm.sq.hat > 0){
    temp = sqrt(diag(diag(Sigma.hat))) %*% (diag(p) - B.hat %*% matrix(c(rep(0,px),A.mat %*% gamma.hat),nrow = 1) / gamma.norm.sq.hat)
  }else{
    temp = sqrt(diag(diag(Sigma.hat)))
  }
  
  SS.hat = temp %*% Omega.hat %*% Sigmae %*% Omega.hat %*%  t(temp)
  SS.Z = SS.hat[(px+1):p,(px+1):p]
  
  set.seed(seed) 
  

  boostrap_sample = MASS::mvrnorm(n = MM, rep(0,pz),SS.Z)
  
  
  sup_sample = rowMaxs(abs(boostrap_sample))
  # I_sup = rowMean(sup_sample)
  
  cv_sup = sort(sup_sample)[floor(MM * (1-sig.level))]
  sup_test = max(abs(sqrt(n)*( sqrt(A.mat) %*% pi.hat.db))) > cv_sup 
  p_val_sup = 1 - ecdf(sort(sup_sample))(max(abs(sqrt(n)*( sqrt(A.mat) %*% pi.hat.db))))
  
  u3 = Omega.hat %*% matrix(c(rep(0,px),A.mat %*% pi.hat),ncol = 1)
  pi.norm.sq.hat <- as.numeric(t(pi.hat) %*% A.mat %*% pi.hat) + 2/n * sum((W %*% u3) * (zeta.hat - W %*% C.hat))
  
  sup_test_power_enhanced = max( c(max(abs(sqrt(n)*( sqrt(A.mat) %*% pi.hat.db))) , sqrt(n)*log(p)*pi.norm.sq.hat )) > cv_sup
  p_val_power_enhanced = 1 - ecdf(sort(sup_sample))( max( c(max(abs(sqrt(n)*( sqrt(A.mat) %*% pi.hat.db))) , sqrt(n)*log(p)*pi.norm.sq.hat )) )
  
  if (beta.output){
    return(list(M_test = sup_test, 
                PM_est = sup_test_power_enhanced, 
                sig.level = sig.level,  
                p_val_M = p_val_sup,
                p_val_PM =  p_val_power_enhanced,
                beta.hat = beta.hat,
                CI = c(CI.low,CI.up)))
  }else{
    return(list(M_test = sup_test, 
                PM_test = sup_test_power_enhanced, 
                sig.level = sig.level,  
                p_val_M = p_val_sup,
                p_val_PM =  p_val_power_enhanced))
  }
  
  
}

source("L1Penalty_functions.R")
