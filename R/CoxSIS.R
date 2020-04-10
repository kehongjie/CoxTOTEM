##' Sure Independence Screening (SIS) for Cox Models with multiple studies.
##' 
##' This function implements the two-step aggregation screening method to Cox 
##' models with multiple studies. It is the first stage (i.e. screening stage)
##' of Cox-TOTEM algorithm.
##' 
##' @param data.list: A list in which every element represents a study. 
##' Within each element (study), \code{time} is the follow up time (n * 1), 
##' \code{event} is the status indicator (n * 1) with 0=alive and 1=dead, and
##' \code{X} is the covariates matrix of dimensions n * p.
##' @param alpha1: Tuning parameter, as \eqn{\alpha_1} in the paper. The Normal
##' critical qunatile when screening within each study. The default is 1e-4.
##' @param alpha2: Tuning parameter, as \eqn{\alpha_2} in the paper. The 
##' Chi-square critical quantile when screening aggregation of studies. The 
##' default is 0.05.
##' 
##' @return A list in which \code{stats} is the vector of the screening 
##' statistics for all the covariates, and \code{index} is the vector of 
##' indices of potential nonzero coefficients.
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
##' }

CoxSIS <- function(data.list, alpha1=1e-4, alpha2=0.05) {
  
  coxz.mat <- multiCoxZ(data.list)
  
  ## step 1. within each study 
  zcut <- qnorm(1-alpha1)
  coxz.mat2 <- coxz.mat  
  coxz.mat2[abs(coxz.mat) >zcut] <- 0
  
  ## step 2. aggregate 
  chi2 <- apply(coxz.mat2,1,function(x) sum((x[x!=0])^2))
  df <- apply(coxz.mat2,1,function(x) sum(x!=0) )
  
  allstrong.index <- which(df==0)
  left.index <- setdiff(1:p,allstrong.index)
  chi2p <- rep(0,length(chi2))
  chi2p[left.index] <-  1-pchisq(chi2[left.index],df[left.index])   
  
  final.index <- union(allstrong.index, which(chi2p < alpha2))
  
  return(list(stats=chi2p,index=final.index))
}


multiCoxZ <- function(data.list){
  p <- ncol(data.list[[1]][["X"]])
  K <- length(data.list)
  coxz.mat <-  matrix(NA,p,K)
  for(k in 1:K){
    dat <- data.list[[k]]
    X <- dat[["X"]]
    time <- c(dat[["time"]])
    event <- c(dat[["event"]])
    
    coxz.mat[,k] <- sapply(1:p, function(j){ 
      datj <- data.frame(time=time,event=event,x=X[,j])
      return(margCoxZ(datj))
    })
  }	
  return(coxz.mat)	
}

margCoxZ <- function(data){ 
  ## data has columns: time, event and x
  c1 <- coxph(Surv(time, event) ~ x , data=data) 
  return(summary(c1)$coef[,"z"])
}




