##' Fit a multiple-study Cox model with group lasso regularization.
##' (or?) Group lasso for Cox model with multiple studies.
##' 
##' This function adopts a group lasso penalty in the partial log-likelihood
##' for Cox model with multiple studies. It is used after the screening by 
##' function \code{CoxSIS} to select the final set of features with
##' nonzero coefficients. The optimization problem is solved by an ADMM
##' algorithm. It is the second stage (i.e. regularization stage) of 
##' Cox-TOTEM algorithm.
##' 
##' @param data.list: A list in which every element represents a study. 
##' Within each element (study), \code{time} is the follow up time (n * 1), 
##' \code{event} is the status indicator (n * 1) with 0=alive and 1=dead, and
##' \code{X} is the covariates matrix of dimensions n * p.
##' @param index: The indices of potential nonzero coefficients selected by
##' screening stage. Typically the "index" output form the function \code{CoxSIS}.  
##' @param lambda: The tuning parameter for regularization. It should be a 
##' single value. 
##' @param rho: TBD.(?) 
##' @param iter: Maximum iterations. Default is 100.
##' @param tol: Convergence threshold for ADMM algorithm. Default is \code{1E-4}.

##' @return A matrix of estimated coefficients, in which each row is a predictor
##' selected by the screening stage and each column is a study.
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
##' res.CoxSIS <- CoxSIS(data.list=data.list, alpha1=1e-4, alpha2=0.05)
##' res.CoxSIS$index
##' 
##' ## group lasso
##' sis.ind <- res.CoxSIS$index
##' coef.final <- CoxGrplasso(data.list,sis.ind,lambda=30,rho=50,iter=100,tol=1e-4)
##' head(coef.final)
##' }


CoxGrplasso <- function(data.list, index, lambda, rho=50, iter=100, tol=1e-4) {
  p <- length(index)
  K <- length(data.list)    
  dat.list <- NULL
  for(k in 1:K){
    dat <- data.list[[k]]
    X <- dat[["X"]]
    X <- X[,index]
    time <- c(dat[["time"]])
    event <- c(dat[["event"]])
    ind <- order(time)        
    datk <- list(time=time[ind],event=event[ind],x=X[ind,])
    dat.list[[k]] <- datk
  }
  ## initialization
  Theta <- matrix(0,nrow=p,ncol=K)
  z.mat <- matrix(0,nrow=p,ncol=K)
  u.mat <- matrix(0,nrow=p,ncol=K)
  ## iterations
  for(i in 1:iter){
    ## update on Theta
    for(k in 1:K){
      Theta[,k] <- updateTheta(dat=dat.list[[k]],theta=Theta[,k],z=z.mat[,k],
                               u=u.mat[,k],rho=rho,tol=tol)
    }
    ## update on Z
    z.mat.old <- z.mat
    for(j in 1:p){
      z.mat[j,] <- updateZ(x=Theta[j,],u=u.mat[j,],lambda=lambda,rho=rho)
    }
    ## update on U
    u.mat <- u.mat+Theta-z.mat
    ## convergence
    if(sum((Theta-z.mat)^2)<tol^2 & sum((z.mat-z.mat.old)^2)<tol^2){
      break
    }
  }
  return(Theta)
}

updateTheta <- function(dat,theta,z,u,rho,iter=100,tol=1e-6){
  n <- dim(dat$x)[1]
  p <- dim(dat$x)[2]
  for(j in 1:iter){
    expXtheta <- as.vector(exp(dat$x%*%theta))
    XexpXtheta <- dat$x*expXtheta
    M0 <- rev(cumsum(rev(expXtheta)))
    M0inv <- ifelse(M0==0,0,1/M0)
    M1 <- NULL
    for(i in 1:p){        
      M1 <- cbind(M1,rev(cumsum(rev(XexpXtheta[,i]))))
    }
    M1M1M0inv2 <- array(0,dim=c(p,p,n)) ## row row^T in M1 /M0^2
    X2 <- array(0,dim=c(p,p,n)) ## XX^T  exp(X^T theta)
    M2 <- array(0,dim=c(p,p,(n+1))) ## sum XX^T exp(X^T theta)
    M2M0inv <- array(0,dim=c(p,p,n)) ## sum XX^T exp(X^T theta)/M0
    for(i in 1:n){
      X2[,,i] <- dat$x[i,]%*%t(dat$x[i,])*expXtheta[i]
      M1M1M0inv2[,,i] <- M1[i,]%*%t(M1[i,])*M0inv[i]^2*dat$event[i]
    }
    for(i in n:1){
      M2[,,i] <- M2[,,(i+1)]+X2[,,i]
      M2M0inv[,,i] <- M2[,,i]*M0inv[i]*dat$event[i]
    }
    ## drv: ell dot with penatly, ddrv: Hessian    
    drv <- -colSums(dat$event*(dat$x-M1*M0inv))+rho*(theta-z+u)
    ddrv  <- rowSums(M2M0inv-M1M1M0inv2,dims=2)+rho*diag(p)
    theta.new <- theta-solve(ddrv)%*%drv
    if(sum(abs(theta.new-theta))<tol){
      break
    }else{
      theta <- theta.new
    }
  }
  return(theta.new)
}

updateZ <- function(x,u,lambda,rho){
  c <- 1-lambda/rho/sqrt(sum((x+u)^2))
  return(ifelse(c>0,c*(x+u),0*x))   
}
