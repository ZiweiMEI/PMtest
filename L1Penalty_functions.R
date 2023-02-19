
fastclime.selector <- function (lambdamtx, icovlist, lambda) 
{
  gcinfo(FALSE)
  d <- dim(icovlist[[1]])[2]
  maxnlambda <- dim(lambdamtx)[1]
  icov <- matrix(0, d, d)
  adaj <- matrix(0, d, d)
  seq <- rep(0, d)
  threshold <- 1e-05
  status <- 0
  for (i in 1:d) {
    temp_lambda <- which(lambdamtx[, i] > lambda)
    seq[i] <- length(temp_lambda)
    if ((seq[i] + 1) > maxnlambda) {
      status <- 1
      icov[, i] <- icovlist[[seq[i]]][, i]
    }
    else {
      icov[, i] <- icovlist[[seq[i] + 1]][, i]
    }
  }
  icov <- icov * (abs(icov) <= abs(t(icov))) + t(icov) * (abs(icov) > 
                                                            abs(t(icov)))
  tmpicov <- icov
  diag(tmpicov) <- 0
  adaj <- Matrix(tmpicov > threshold, sparse = TRUE) * 1
  sparsity <- (sum(adaj))/(d^2 - d)
  # if (status == 1) {
  #   cat("Some columns do not reach the required lambda!\n You may want to increase lambda.min or use a larger nlambda. \n")
  # }
  rm(temp_lambda, seq, d, threshold)
  gc()
  result <- list(icov = icov, adaj = adaj, sparsity = sparsity)
  class(result) = "fastclime.selector"
  return(result)
}


get_lasso_lambda_max <- function(x, y, scalex = FALSE, nlambda = 100, lambda_min_ratio = 0.0001){
  
  n <- nrow(x)
  p <- ncol(x)
  y <- as.numeric(y)
  
  sd_n <- function(x){
    sqrt(sum((x - mean(x))^2)/length(x))
  }
  
  
  if (scalex) {
    sx <- scale(x, scale = apply(x, 2, sd_n))
    sx <- as.matrix(sx, ncol = p, nrow = n)
    
    lambda_max <- max(abs(colSums(sx * y))) / n
  } else {
    lambda_max <- max(abs(colSums(x * y))) / n
  }
  
  return(lambda_max)
  
}

get_lambda_seq <- function(lambda_max, lambda_min_ratio = 0.001, nlambda = 100){
  
  lambda_min <- lambda_min_ratio * lambda_max
  lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
  
  return(lambda_seq)
}

Lasso.bic <- function(x,y,standardize,intercept = T,weight = rep(1,ncol(x))){
  
  
  lambda_max <- get_lasso_lambda_max(x,y,standardize)
  lambda_seq <- get_lambda_seq(lambda_max)
  
  n = dim(x)[1]
  
  out <- glmnet(x ,y, standardize = standardize, lambda = lambda_seq, penalty.factor = weight)
  df <- out$df + 1
  
  coef.Las <- coef(out,s = out$lambda)
  y.hat <-  cbind(1,x) %*% coef.Las 
  mse <- colMeans((y.hat - matrix(y,n,length(out$lambda)))^2)
  BIC <- n * mse + log(n) * df
  # print(which.min(BIC))
  lambda.bic <- out$lambda[which.min(BIC)]
  out <- glmnet(x ,y, standardize = standardize, lambda = lambda.bic, penalty.factor = weight)
  coef.bic <- coef(out,s = out$lambda.bic)
  return(coef.bic)
}


foldid_vec = function(TT, k) {
  
  seq.interval = split(1:TT, ceiling(seq_along(1:TT)/(TT/k)))
  
  id = rep(0, TT)
  for (j in 1:k) {
    id[seq.interval[[j]]] = j
  }
  
  return(id)
}



Lasso <- function(x,y,standardize = T,intercept = T, tuning = "CV", weight = rep(1,ncol(x)), block.CV = TRUE, k = 10, seed = 2023){
  
  set.seed(seed)
  
  if ( is.null(tuning) ){
    if (block.CV){
      foldid = foldid_vec(nrow(x),k)
      out.Lasso <- cv.glmnet(x,y,standardize = standardize,intercept = intercept, penalty.factor = weight,
                               foldid = foldid)
      
    }else{
      out.Lasso <- cv.glmnet(x,y,standardize = standardize,intercept = intercept, penalty.factor = weight)
    }
    
    coef <-  as.vector(coef(out.Lasso, s = out.Lasso$lambda.min)) 
  }else if (tuning == "CV"){
    if (block.CV){
      foldid = foldid_vec(nrow(x),k)
      out.Lasso <- cv.glmnet(x,y,standardize = standardize,intercept = intercept, penalty.factor = weight,
                             foldid = foldid)
    }else{
      out.Lasso <- cv.glmnet(x,y,standardize = standardize,intercept = intercept, penalty.factor = weight)
    }
    coef <-  as.vector(coef(out.Lasso, s = out.Lasso$lambda.1se)) 
  }else if (tuning == "BIC"){
    coef <- Lasso.bic(x,y,standardize,intercept, weight = weight)
  }
  
  if (!intercept){
    coef <- coef[-1]
  }
  
  coef = matrix(coef)
  
  return(coef)

}