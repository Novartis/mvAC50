
#################################################################################
# Copyright 2019 Novartis Institutes for BioMedical Research Inc.               #
# Licensed under the Apache License, Version 2.0 (the "License");               #
# you may not use this file except in compliance with the License.              #
# You may obtain a copy of the License at                                       #
#                                                                               #
# http://www.apache.org/licenses/LICENSE-2.0                                    #
#                                                                               #
# Unless required by applicable law or agreed to in writing, software           #
# distributed under the License is distributed on an "AS IS" BASIS,             #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.      #
# See the License for the specific language governing permissions and           #
# limitations under the License.                                                #
#################################################################################



#********************************************************
# Selected Dose Response Curve Analysis (Core) Functions
#********************************************************
# Author:
# Hanspeter Gubler
# Novartis Institutes for BioMedical Research
# NIBR Informatics / IS Delta
# WSJ-310.5.17
# Basel, Switzerland
#********************************************************
#
# File    : DRC Analysis Core Functions
# Version : V1.2.2
# Date    : June 8, 2019
#
#********************************************************
# Purpose : 
# 
# Provide a set of functions to support fitting of different model curves to dose-response data:
#
# - 4 parameter logistic (4PL) / Hill curve fitting (sigmoidal)
# - nonparametric local regression curve fitting to characterize non-sigmoidal curves and to provide optimal starting 
#   parameters for 4PL fitting
# - constant fitting for the characterizaion of data without / weak dose-dependence of responses --> determination of
#   average response
#
#********************************************************
# History :
#
# V1.2.2 June 8, 2019:
#        Code cleanup for publication in github
# V1.2.1 July 24, 2017, August 8, 2017
#        changes by Steffen Renner
#        bug fix: x replaced with log.x in nonparametric.fit.results function (determination of a fallback value for ac50)
#        4PL fit function extended with constraints for A0 and Ainf, all constraints now passed as arguments to
#        param.4PL.fit.results() function
# V1.2.0 July 19, 2017:
#        Functions assembled from various other scripts and streamlined to unify the return data structure for nonparametric, 
#        4PL/Hill-curve and constant fits.
#
#********************************************************
#
# Notes (Restrictions) for all V1.2.x script functions:
#
# - High level multi-model analysis logic is principally set up for 'positive-going' responses. Some 
#   criteria and comparison operations are not  equipped or not tested for descending/negative-going curves.
#
# - Standard error estimates of fitted parameters are not derived in this version of the core function set.
# 
#   These can be easily added for 4PL and constant fits. Approximate pointwise confidence bands for the fitted 
#   parametric (nonlinear) 4PL can then e.g. be based on error propagation (delta rule) aproach and is trivial to do for the
#   constant model. 4PL confidence bands must be viewed with caution when parameter correlations are high and/or parameter
#   standard errors are large. The confidence band  for the nonparametric locfit model can be generated via suitable 
#   'predict(model, x values, se.fit=TRUE)' calls. 
#
# - Fits of bell-shaped curves are presently only supported when using the nonparametric fitting approach. Parametric 
#   bell-curve models (using additive or multiplicative terms) can be added in future, if such a need arises.
#
#********************************************************
# operating system and server properties

Sys.getlocale()
.Platform

#********************************************************
# Required packages

# nonparametric (local regression) fitting
library(locfit)

# robust linear regression, rlm()
library(MASS)

# robust nonlinear regression, nlrob()
library(robustbase)


#********************************************************
nonparametric.fit <- function(x, y, deg=2, alpha=0.4)
{
  f <- f1 <- f2 <- list(fit=NULL, model=NULL)
  
  f  <- locfit(y ~ x, deg=deg, alpha=alpha, lfproc=locfit.raw)
  
  if (deg > 0) {  
    f1 <- first.deriv(x, y, alpha=alpha)
	if (deg > 1) {
      f2 <- second.deriv(x, y, alpha=alpha)
	} 
  } 
  
  list(fit=fitted(f), model=f, 
       fder=f1$fit, fder.model=f1$model,
       sder=f2$fit, sder.model=f2$model)
  
} # nonparametric.fit

#********************************************************
first.deriv <- function(x, y, deg=2, alpha=0.4)
{
  f <- locfit(y ~ x, deg=deg, alpha=alpha, deriv=1, lfproc=locfit.raw) 
  list(fit=fitted(f), model=f)
  
} # first.deriv

#********************************************************
second.deriv <- function(x, y, deg=2, alpha=0.4)
{
  f <- locfit(y ~ x, deg=deg, alpha=alpha, deriv=c(1,1), lfproc=locfit.raw) 
  list(fit=fitted(f),model=f)
  
} # second.deriv

#********************************************************
values <- function(model, x)
{
  if ((length(x) == 1)&&(is.na(x))) return(NA)
  predict(model, newdata=x)
}

#********************************************************
# when model is a first derivative model (fder.model) 
# slopes(model, x)

slopes <- values

#********************************************************
# single trapezoid area (pairs of x,y)

trap.area <- function(x, y)
{
  if ((length(x) != 2)||(length(y) !=2)) return(NA)
  
  0.5*(y[2] + y[1])*(x[2] - x[1])
  
} # trap.area

#********************************************************
# trapezoidal integration over series of x, y

trap.integ <- function(x,y) 
{
  if ((length(x) <= 1)||(length(y) <= 1)) return(NA)  
  integ <- 0
  for (i in 2:length(x)) {
    integ <- integ + trap.area(c(x[i-1], x[i]), c(y[i-1], y[i]))    
  }
  
  integ
  
} # trap.integ 

#********************************************************
# 4PL, 4-Paramater Logistic (Hill-) Curve Function

hill.curve.logscale <- function(a0, ainf, logac50, hillslope, logx)
{
  ainf + (a0-ainf)/(1.0 + 10^(hillslope*(logx - logac50))) 
}

#********************************************************
# Hill Curve Inversion
#********************************************************
# this uses ac50, not log(ac50) !

inverted.hill.curve <- function(a0, ainf, ac50, hillslope, y)
{
  ac50*((a0-y)/(y-ainf))^(1/hillslope)
}

#********************************************************
# Correlation Matrix from Covariance Matrix
#********************************************************
corr <- function(covar)
{ 
  d <- dim(covar)
  if (d[1] != d[2]) return(NULL)
  
  nd <- d[1]   
  corr <- matrix(NA, nd, nd)
  
  for (i in 1:nd) {
    for (j in 1:nd) {
      corr[i,j] <- covar[i,j]/sqrt(covar[i,i]*covar[j,j])
    }	  
  }
  corr
}

#********************************************************
# Numerical Root Finding, find *all* x where f(x)=0
# where x and y are given as arrays
#********************************************************
find.roots <- function(x, y, type)
{
  # type 1: pos -> neg transition root locations
  # type 2: neg -> pos transition root locations
  # type 3: both transition directions at once
  
  n          <- length(y)
  
  isect      <- rep(NA, n)
  isect.type <- rep(NA, n)
  
  num.isect <- 0
  
  for (i in 1:(n-1)) {
    
    y1 <- y[i]
    y2 <- y[i+1]
    
    if ((type == 1)||(type== 3)) {
      if ((y2 < 0)&&(y1 >= 0)) {
        num.isect <- num.isect+1
        x1 <- x[i]
        x2 <- x[i+1]
        isect[num.isect] <- x1 - (x2-x1)/(y2-y1)*y1
        isect.type[num.isect] <- 1
      }
    } 
    
    if ((type == 2)||(type == 3)) {
      if ((y2 > 0)&&(y1 <= 0)) {
        num.isect <- num.isect+1
        x1 <- x[i]
        x2 <- x[i+1]
        isect[num.isect] <- x1 - (x2-x1)/(y2-y1)*y1
        isect.type[num.isect] <- 2           
      }     
    }
    
  } # for
  
  # truncate arrays to appropriate length
  if (num.isect > 0) {
    isect      <- isect[1:num.isect]
    isect.type <- isect.type[1:num.isect]
  } else {
    isect      <- isect[1]
    isect.type <- isect.type[1]
  }  
  
  list(num=num.isect, isect.type=isect.type, isect=isect)
  
} # find.roots

#********************************************************
# parametric (4PL) fit and result assembly
# need to obtain initial parameter estimates ('params' argument) from a prior nonparametric fit
#********************************************************

param.4PL.fit.results <- function(logx, y, params, method="robust", abs.acx.level=50, constr.lower, constr.upper)
  {
  # signal amplitude threshold to check for a0/ainf differences
  eps <- 1E-6
  
  fit <- NULL
  
  n.data   <- length(logx)
  rweights <- rep(1, n.data)
  
  # initial parameter estimates are set from nonparametric estimates and passed into 'params' array argument by the caller
  
  if (length(params) !=4) {
    
    return(list(rse=NA, r.squared=NA, a0=NA, ainf=NA, ac50=NA, hill=NA, abs.acx=NA, abs.acx.level=abs.acx.level, 
           convergence=-1, message="inital parameter estimates were not properly set", rweights=NA,
           fit.type="param.4PL", fit=NULL)) 
    
  } else {
    
    a0.init   <- params[1]
    ainf.init <- params[2]
    if (abs(a0.init - ainf.init) < eps) ainf.init <- a0.init + (max(y) - min(y)) 
    ac50.init <- params[3]
    if ((ac50.init > max(10^logx))||(ac50.init < min(10^logx)))  ac50.init <- 10^mean(logx)
    hill.init <- params[4]
    if ((hill.init > 10)||(hill.init < 0.1)) hill.init <- 1    
  }  

  if (method == "robust") {

    require(robustbase)
    
    try(fit <- nlrob(y ~ hill.curve.logscale(a0, ainf, logac50, hill, logx), 
                   data=data.frame(logx, y),   
                   method="M",
                   algorithm="port",
                   start = c(a0=a0.init, ainf=ainf.init, logac50=log10(ac50.init), hill=hill.init),
                   lower = constr.lower,
                   upper = constr.upper,
                   tol=1E-5,	  # default: 1E-6,
                   maxit=100,	  # default: 20
                   doCov=TRUE,                   
                   trace = FALSE,
                   control=list(warnOnly = TRUE)
                  )
        )

    fit.type <- "param.4PL/nlrob"            
        
    # switch the method if robust fit failed and try again (below)
    if (is.null(fit)) method <- "nls"
  }
  
  if (method != "robust") {
    
    try(fit <- nls(y ~ hill.curve.logscale(a0, ainf, logac50, hill, logx), 
                    algorithm="port",
                    start = list(a0=a0.init, ainf=ainf.init, logac50=log10(ac50.init), hill=hill.init),
                    lower = list(constr.lower[1], constr.lower[2], constr.lower[3], constr.lower[4]),
                    upper = list(constr.upper[1],  constr.upper[2],  constr.upper[3],  constr.upper[4]),
                    trace = FALSE,
                    control=list(warnOnly = TRUE)
                  )
        )
    
     fit.type <- "param.4PL/nls"
  }  
  
  if (!is.null(fit)) {  	# do this only if the fit succeeded	...
    # summary(fit)

    # residual standard error 
    rse <- sqrt(deviance(fit)/(n.data-4))
    
    # fit diagnostics
    
    if (fit.type == "param.4PL/nls") {
      convergence <- fit$convergence
      message     <- fit$message
      
    } else {
      convergence <- fit$status
      message     <- fit$status
      rweights    <- fit$rweights
    }
    
    pars <- as.numeric(coef(fit))
    attr(pars, "names") <- rep(NULL, 4)

    a0   <- pars[1]
    ainf <- pars[2]        
    ac50 <- 10^pars[3]
    hill <- pars[4]

    # this can be handled the same way also for non-monotone functions
    log.x       <- seq(min(logx), max(logx), length.out=101)        
    fitted.line <- predict(fit, list(logx=log.x))
    
    amax		<- max(fitted.line)
    log.c.amax	<- log.x[which.max(fitted.line)]
    c.amax		<- 10^log.c.amax
    
    amin		<- min(fitted.line)
    log.c.amin	<- log.x[which.min(fitted.line)]
    c.amin		<- 10^log.c.amin
    
    # here we determine abs.acx numerically instead of through function inversion
    
    # decreasing or increasing curve ?
    if (ainf > a0) {
      # ascending curve
      tmp <- find.roots(log.x, fitted.line - abs.acx.level, 2)
    } else {
      # descending curve
      tmp <- find.roots(log.x, fitted.line - abs.acx.level, 1)
    }  
    
    if (tmp$num > 0) {
      abs.acx <- 10^tmp$isect[1]
    } else {
      abs.acx <- NA_real_
    }            

    # residual sum of squares
    ss.res <- sum((y - predict(fit))^2)
    ss.res	

    # constant for deriving total ss
    const.fit <- mean(y)
    const.fit

    ss.tot <- sum((y - const.fit)^2)
    ss.tot
    
    # R^2, fraction of explained variance
    r.squared <- 1 - ss.res/ss.tot
    r.squared	
      
  } else {
    # fit did not succeed ...    
    # fill in default results
    
    rse  <- r.squared <- NA
    a0   <- ainf <- ac50 <- abs.acx <- hill <- NA
    amax <- c.amax <- amin <- c.amin <- NA
    convergence  <- -1
    message      <- "fit failed"
    fit          <- NULL
    
  }
  
  res <- list(rse=rse, r.squared=r.squared, a0=a0, ainf=ainf, ac50=ac50, hill=hill, abs.acx=abs.acx, abs.acx.level=abs.acx.level,
              amax=amax, c.amax=c.amax, amin=amin, c.amin=c.amin,              
              convergence=convergence, message=message, rweights=rweights, fit.type=fit.type, fit=fit)
  res
  
} # param.4PL.fit.results


#********************************************************
# nonparametric fit and result assembly
#********************************************************

nonparametric.fit.results <- function(logx, y, abs.acx.level=50)
{
  # used as a signal amplitide threshold for hill-slope estimation
  eps <- 1E-6
  
  require(locfit)
  
  n.conc <- length(unique(logx))
  n.data <- length(logx)
  # average # of replicates at each conc  
  n.rep  <- n.data / n.conc      
  
  # heuristic choice for 'optimal' alpha in [0.2, 1]
  alpha  <- min(max(3.25/n.conc, 0.2), 1)

  locreg <- nonparametric.fit(logx, y, deg=1, alpha=alpha)
  
  # residual sum of squares
  ss.res <- sum((y - values(locreg$model, logx))^2)
  ss.res
  
  # residual stanadard error, degrees of freedom not considered
  # rse  <- sqrt(ss.res / n.data)
  
  # residual scale with consideration of proper number of degrees of freedom
  rse    <- as.numeric(sqrt(locreg$model$dp["rv"]))
  
  const.fit <- mean(y)
  const.fit
  
  ss.tot <- sum((y - const.fit)^2)

  # R^2, fraction of explained variance using spline fit
  r.squared <- 1 - ss.res/ss.tot

  # rough curve characterization ...
  a0      <- values(locreg$model, min(logx))
  ainf    <- values(locreg$model, max(logx))

  # determine ~ ac50, log.x (finer grid)
  
  log.x       <- seq(min(logx), max(logx), length.out=101)          
  fitted.line <- values(locreg$model, log.x)

  id.max  <- which.max(values(locreg$model, log.x))
  logxmax <- log.x[id.max]
  amax    <- values(locreg$model, logxmax)
  
  id.min  <- which.min(values(locreg$model, log.x))
  logxmin <- log.x[id.min]
  amin    <- values(locreg$model, logxmin)    
  
  
  # determine principal directionality of curve
  df <- data.frame(logx, y)
  m  <- aggregate(df$y, list(df$logx), mean)
  
  AUC <- trap.integ(m[,1], m[,2])
  
  if (AUC > 0) {
    if (amax > ainf) {
      ref.level <- amax
    } else {
      ref.level <- ainf
    }  
  } else {
    if (amin < ainf) {
      ref.level <- amin
    } else {
      ref.level <- ainf
    }
  }

  ac50.level     <- 0.5*(a0 + ref.level)
  
  # relative ac50 estimate
  if (ref.level > a0) {
    # ascending curve
    tmp <- find.roots(log.x, fitted.line - ac50.level, 2)
  } else {
    # descending curve
    tmp <- find.roots(log.x, fitted.line - ac50.level, 1)
  }  
  
  if (tmp$num > 0) {
    ac50 <- 10^tmp$isect[1]
  } else {
    # corrected in V1.2.1
    ac50 <- 10^(mean(log.x))
  }           
  
  # absolute ac50 estimate
  if (ref.level > a0) {
    # ascending curve
    tmp <- find.roots(log.x, fitted.line - abs.acx.level, 2)
  } else {
    # descending curve
    tmp <- find.roots(log.x, fitted.line - abs.acx.level, 1)
  }  
  
  if (tmp$num > 0) {
    abs.acx <- 10^tmp$isect[1]
  } else {
    abs.acx <- NA_real_
  }             
  
  # determine approximate hill slope equivalent
  
  fder.estimate <- slopes(locreg$fder.model, log10(ac50))
  
  delta.1 <- abs(a0 - ainf)
  delta.2 <- abs(a0 - amax)
  delta   <- max(delta.1, delta.2)
  
  if (delta > eps ) { 
    hill <- 4*fder.estimate/(log(10)*delta)
  } else {
    hill <- NA_real_
  }  
  
  res <- list(rse=rse, r.squared=r.squared, 
              a0=a0, ainf=ainf, ac50=ac50, hill=hill, abs.acx=abs.acx, abs.acx.level=abs.acx.level,
              amax=amax, c.amax=10^logxmax, amin=amin, c.amin=10^logxmin,
              convergence="", message="", alpha=alpha, fit.type="nonpar", fit=locreg$model)
  res  
  
} # nonparametric.fit.results


#********************************************************
# constant fit and result assembly
# fill in 'standardized' data structure to make equivalent to parametric result structure
#********************************************************

constant.fit.results <- function(logx, y, method="robust", abs.acx.level=NA_real_ )
{
  
  n.data <- length(logx)

  if (method != "robust") {
    
    fit <- lm(y ~ 1)
    pred <- predict(fit)
    message  <- "mean"
    fit.type <- "lm"
    
  } else if (method == "robust") {

    require(MASS)
    
    fit <- rlm(y ~ 1)
    pred <- as.numeric(predict(fit))
    message  <- "robust (biweight) mean"
    fit.type <- "rlm/robust"
  }  

  # residual sum of squares
  ss.res <- sum((y - pred)^2)
  rse    <- sqrt(ss.res / (n.data-1) )

  # R^2, fraction of explained variance
  r.squared <- 0

  res <- list(rse=rse, r.squared=r.squared,
              a0=pred[1], ainf=pred[n.data], ac50=NA_real_, hill=0, abs.acx=NA_real_, abs.acx.level=NA_real_,
              amax=pred[1], c.amax=NA_real_, amin=pred[1], c.amin=NA_real_,
              convergence="", message=message, fit.type=fit.type, fit=fit)
  res  
  
} # constant.fit.results


#********************************************************
# END Dose Response Curve Analysis Core Functions
#********************************************************
	