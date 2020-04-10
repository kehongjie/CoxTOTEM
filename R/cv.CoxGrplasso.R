##' Cross-validation for CoxGrplasso.
##'  
##' This function does M-fold cross-validation for \code{CoxGrplasso}.
##' 
##' @param data.list: A list in which every element represents a study. 
##' Within each element (study), \code{time} is the follow up time (n * 1), 
##' \code{event} is the status indicator (n * 1) with 0=alive and 1=dead, and
##' \code{X} is the covariates matrix of dimensions n * p.
##' @param index: The indices of potential nonzero coefficients selected by
##' screening stage. Typically the "index" output form the function \code{CoxSIS}.  
##' @param rho: TBD.(?) 
##' @param iter: Maximum iterations. Default is 100.
##' @param tol: Convergence threshold for ADMM algorithm. Default is \code{1E-4}.
##' @param tlmbd: A vector of user-supplied optional lambda sequence.
##' @param M: Number of folds in cross-validation. Default is 5.
##' 
##' @return A list with two elements: \code{rate} is the mean rate for each lambda
##' value and \code{lambda.final} is the optimal lambda value selected by 
##' cross-validation.
##' @export
##' @examples
##' \dontrun{
##' set.seed(1234)
##' library(MASS)
##' library(survival)
##' 
##' n <- 50 ## sample size
##' p <- 200 ## number of covariates
##' rho <- 0.5
##' K <- 5 ## number of studies
##' ss <- 2 ## number of signal/true predictors
##' lambda0 <- 1 ## baseline hazard
##' rate <- 0.2 ## parameter for Exponential distribution
##' mu <- 1 ## signal size
##' 
##' ## set up the coefficients
##' true.ind <- sample(1:p,size=ss) ## signal index
##' noise.ind <- (1:p)[-true.ind]
##' beta.mat <- matrix(0,nrow=p,ncol=K)
##' for(jj in 1:length(true.ind)){
##'   beta.mat[true.ind[jj],] <- mu
##' }
##' 
##' ## simulate data for each study
##' data.list <- vector("list",K)
##' for(k in 1:K){
##'   beta <- beta.mat[,k]
##'   sigma <- toeplitz(rho^c(0,1:(p-1)))
##'   X <- mvrnorm(n, rep(0,p), sigma)
##'   U <- runif(n,0,1)
##'   C <- rexp(n,rate=rate)
##'   time <- -log(U) / (lambda0*exp(X%*%beta)) 
##'   Y <- pmin(time,C)
##'   D <- ifelse(time<C,1,0)
##'   data.list[[k]] <- list(time=Y,event=D,X=X)
##' }
##' 
##' ## screening
##' res.SIS <- CoxSIS(data.list=data.list, alpha1=1e-4, alpha2=0.05)
##' res.CoxSIS$index
##' 
##' ## group lasso cross-validation
##' sis.ind <- res.CoxSIS$index
##' tlmbd <- c(30,35,38,40) ## optional lambda sequence
##' res.cv <- cv.CoxGrplasso(data.list, index=sis.ind, rho=50, iter=100, tol=1e-4, tlmbd, M=5)
##' print(res.cv$rate)
##' print(res.cv$lambda.final)
##' }


cv.CoxGrplasso <- function(data.list, index, rho=50, iter=100, tol=1e-4, tlmbd, M=5) {
  # GLASSO
  coef.lasso <- list()
  sens.lasso <- vector()
  spec.lasso <- vector()
  fnllss.ind <- list()

  ## split the data (M-fold)
  itest <- list()
  itrain <- list()
  for (i in 1:length(data.list)) {
    dat.evt <- list(event = data.list[[i]]$event[data.list[[i]]$event== 1],
                    time = data.list[[i]]$time[data.list[[i]]$event== 1],
                    X = data.list[[i]]$X[data.list[[i]]$event== 1,])
    dat.csr <- list(event = data.list[[i]]$event[data.list[[i]]$event== 0],
                    time = data.list[[i]]$time[data.list[[i]]$event== 0],
                    X = data.list[[i]]$X[data.list[[i]]$event== 0,])
    cv.ind.evt <- sample(length(dat.evt$event))
    cv.ind.csr <- sample(length(dat.csr$event))
    folds.evt <- cut(seq(1, length(dat.evt$event)), breaks = M, labels = F)
    folds.csr <- cut(seq(1, length(dat.csr$event)), breaks = M, labels = F)
    datcv.evt <- list(time = dat.evt[["time"]][cv.ind.evt], 
                      event = dat.evt[["event"]][cv.ind.evt], X = dat.evt[["X"]][cv.ind.evt, ])
    datcv.csr <- list(time = dat.csr[["time"]][cv.ind.csr], 
                      event = dat.csr[["event"]][cv.ind.csr], X = dat.csr[["X"]][cv.ind.csr, ])
    inter.train <- list()
    inter.test <- list()
    for(j in 1:M){
      tind.evt <- which(folds.evt == j, arr.ind = T)
      tind.csr <- which(folds.csr == j, arr.ind = T)
      test.evt <- list(time = datcv.evt[["time"]][tind.evt],
                       event = datcv.evt[["event"]][tind.evt],X = datcv.evt[["X"]][tind.evt,])
      test.csr <- list(time = datcv.csr[["time"]][tind.csr],
                       event = datcv.csr[["event"]][tind.csr],X = datcv.csr[["X"]][tind.csr,])
      
      train.evt <- list(time = datcv.evt[["time"]][-tind.evt],
                        event = datcv.evt[["event"]][-tind.evt],X = datcv.evt[["X"]][-tind.evt,])
      train.csr <- list(time = datcv.csr[["time"]][-tind.csr],
                        event = datcv.csr[["event"]][-tind.csr],X = datcv.csr[["X"]][-tind.csr,])
      
      inter.test[[j]] <- list(time = c(test.evt[["time"]],test.csr[["time"]]),
                              event = c(test.evt[["event"]],test.csr[["event"]]), 
                              X = rbind(test.evt[["X"]],test.csr[["X"]]))
      inter.train[[j]] <- list(time = c(train.evt[["time"]],train.csr[["time"]]),
                               event = c(train.evt[["event"]],train.csr[["event"]]), 
                               X = rbind(train.evt[["X"]],train.csr[["X"]]))
    }
    itest[[i]] <- inter.test
    itrain[[i]] <- inter.train
  }
  test <- vector("list",length = M)
  train <- vector("list",length = M)
  for(jjj in 1:M){
    for(ii in 1:length(data.list)){
      test[[jjj]][[ii]] <- itest[[ii]][[jjj]]
      train[[jjj]][[ii]] <- itrain[[ii]][[jjj]]
    } 
  }
  
  ## train and test
  rate.cv <- vector()
  coef.glasso <- list()
  glasso.ind <- list()
  for(l in 1:length(tlmbd)){
    ## train
    coef.glasso[[l]] <- lapply(train, FUN = function(x){CoxGrplasso(x, index, tlmbd[l], rho=rho, tol=tol)})
    glasso.ind[[l]] <- lapply(coef.glasso[[l]], FUN = function(x){index[which(rowSums(abs(x[, 1:K]) > tol) == K)]})
    for(s in 1:M){
    }
    ## test
    ratt <- vector()
    for(m in 1:M){
      if(length(glasso.ind[[l]][[m]]) == 0){
        ratt[m] <- 0
      }else{
        condi <- which(rowSums(abs(coef.glasso[[l]][[m]][, 1:K]) > tol) == K)
        ratt[m] <- crossv(test[[m]], glasso.ind[[l]][[m]], coef.glasso[[l]][[m]][condi,]) 
      }
    }
    rate.cv[l] <- mean(ratt,na.rm = T)
  }

  ## lambda chosen
  lambda.final <- tlmbd[which.max(rate.cv)]
  return(list(rate=rate.cv, lambda.final=lambda.final))
}



crossv <- function(data.list,index,theta,t=0){   
  K <- length(data.list)
  t.med <- numeric(K)
  success <- numeric(K)
  if(length(index) == 1){theta = t(theta)}
  for(k in 1:K){        
    dat <- data.list[[k]]
    p <- length(index)
    X <- dat[["X"]]
    X <- as.matrix(X[,index])
    time <- c(dat[["time"]])
    event <- c(dat[["event"]])
    N <- length(event)
    ind <- order(time)        
    dat <- list(time=time[ind],event=event[ind],x=X[ind,])        
    t.med <- median(time)
    ## compute Breslow Estimator
    if(length(index) == 1){
      #theta = t(as.matrix(theta))
      dat$x = as.matrix(dat$x)
    }
    expXtheta <- as.vector(exp(dat$x%*%theta[,k]))
    XexpXtheta <- dat$x*expXtheta
    M0 <- rev(cumsum(rev(expXtheta)))/N
    M0inv <- ifelse(M0==0,0,1/M0)
    ## Breslow estimator (little lambda)
    lambda <- dat$event*M0inv/N
    Lambda <- cumsum(lambda)        
    ## t to check predicton at 
    if(t==0){
      t <- t.med
    }
    ind.t <- which.min(abs(dat$time-t))
    ## Survival function
    if(dat$time[ind.t]-t<0){
      estS.t <- exp(-expXtheta*Lambda[ind.t])
    }else{
      estS.t <- exp(-expXtheta*Lambda[(ind.t-1)])
    }
    ## subset to non-censored subjects
    ind.uc <- (1:N)[dat$event==1]
    m <- length(ind.uc)
    time.uc <- dat$time[ind.uc]
    estS.t <- estS.t[ind.uc]
    ## proportion of correct prediction
    success[k]  <- (sum(time.uc>=t & estS.t>=.5) + sum(time.uc<t & estS.t<.5) )/m
    
  }
  ## average correct prediction rate over studies
  rate <- mean(success)
  return(rate)
}

