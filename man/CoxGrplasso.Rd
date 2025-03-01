% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoxGrplasso.R
\name{CoxGrplasso}
\alias{CoxGrplasso}
\title{Fit a multiple-study Cox model with group lasso regularization.}
\usage{
CoxGrplasso(data.list, lambda, rho = 50, iter = 100, tol = 1e-04)
}
\arguments{
\item{data.list:}{A list in which every element represents a study,
typically after screening. Within each element (study), \code{time} is
the follow up time (\eqn{n \times 1}), \code{event} is the status indicator
(\eqn{n \times 1}) with 0=alive and 1=dead, and \code{X} is the covariates
matrix of dimensions \eqn{n \times p}, where p is the number of covariates
left after screening.}

\item{lambda:}{The penalty parameter in regularization. It is a single
optimal value that needs to be tuned in \code{cv.CoxGrplasso}.}

\item{rho:}{The step size in ADMM algorithm. Default is 50, and change of
\code{rho} is typically not recommended.}

\item{iter:}{Maximum iterations. Default is 100.}

\item{tol:}{Convergence threshold for ADMM algorithm. Default is \code{1E-4}.}
}
\value{
A matrix of estimated coefficients, in which each row is a covariate
and each column is a study.
}
\description{
This function adopts a group lasso penalty to the partial log-likelihood
for Cox model with multiple studies. It is used after the screening by
function \code{CoxSIS} to select the final set of features with
nonzero coefficients. The optimization problem is solved by an ADMM
algorithm. It is the second stage (i.e. regularization stage) of
Cox-TOTEM method.
}
\examples{
\dontrun{
set.seed(1234)
library(MASS)
library(survival)

n <- 50 ## sample size
p <- 200 ## number of covariates
rho <- 0.5
S <- 5 ## number of studies
ss <- 2 ## number of signal/true predictors
lambda0 <- 1 ## baseline hazard
rate <- 0.2 ## parameter for Exponential distribution
mu <- 1 ## signal size

## set up the coefficients
true.ind <- sample(1:p,size=ss) ## signal index
noise.ind <- (1:p)[-true.ind]
beta.mat <- matrix(0,nrow=p,ncol=S)
for(jj in 1:length(true.ind)){
  beta.mat[true.ind[jj],] <- mu
}

## simulate data for each study
data.list <- vector("list",S)
for(s in 1:S){
  beta <- beta.mat[,s]
  sigma <- toeplitz(rho^c(0,1:(p-1)))
  X <- mvrnorm(n, rep(0,p), sigma)
  U <- runif(n,0,1)
  C <- rexp(n,rate=rate)
  time <- -log(U) / (lambda0*exp(X\%*\%beta)) 
  Y <- pmin(time,C)
  D <- ifelse(time<C,1,0)
  data.list[[s]] <- list(time=Y,event=D,X=X)
}

## screening
res.CoxSIS <- CoxSIS(data.list=data.list, alpha1=1e-4, alpha2=0.05)

## remove variables with zero coefficients
index.sis <- res.CoxSIS$index
data.list.sis <- NULL
for(s in 1:S){
 dat <- data.list[[s]]
 X <- dat[["X"]]
 X <- X[,index.sis]
 dat[["X"]] <- X
 data.list.sis[[s]] <- dat
}

## group lasso
coef.final <- CoxGrplasso(data.list.sis,lambda=30,rho=50,iter=100,tol=1e-4)
head(coef.final)
}
}
