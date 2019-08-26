
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


library(tidyverse)

#' Fit AC50s to dose response data
#'
#' Function to calculate best fits for dose response data. Three fits are attempted, constant, 4-parameter logistic regression, 
#' and a nonparametric fit. The best fit is selected based on a heuristic.
#' 
#'
#' @param df dataframe containing concentration and readout columns
#' @param conc_um treatment concentration in uM
#' @param readout_value_pctAC readout value as percentage of an active control
#' @param best_fit by default, only the best fit is reported. If set FALSE, all three fits are reported (the best fit is reported with fit_sel=1 )
#' @param abs_acx_level threshold for absolute AC50s, default is 50 percent activity
#' @param a_diff_thresh if amax-amin < a_diff_thresh, fit is considered constant 
#' @param rsq_threshold r squared threshold for parametric and non-parametric fits to be considered
#' @param constr_lower list of lower limit constraits for parametric fit for a0, ainf, ac50, hill, default c(-50, -50, log10(1e-06), 0.1). Use -Inf for no constraints.
#' @param constr_upper list of upper limit constraits for parametric fit for a0, ainf, ac50, hill, default c( 500,  500,  3,  10). Use Inf for no constraints.
#'
#' @return  data.frame with AC50s and fitted parameters
#' @export
#'
#' @examples
#' see vignettes
#' 
#' 
fit_AC50 <- function(df, conc_um, readout_value_pctAC, best_fit=TRUE, abs_acx_level=50, a_diff_thresh=30, rsq_threshold=0.5,
                    constr_lower=c(-50, -50, log10(1e-06), 0.1), constr_upper=c( 500,  500,  3,  10)){
  
  if (conc_um != "conc_um"){df$conc_um <- df[[conc_um]]}
  if (readout_value_pctAC != "readout_value_pctAC"){df$readout_value_pctAC <- df[[readout_value_pctAC]]}
  
  # filter out potential NAs in readout values
  df <- df %>% filter(!is.na(readout_value_pctAC) & !is.na(conc_um))
  
  abs.acx.level<- abs_acx_level
  a.diff.thresh <-  a_diff_thresh
  rsq.threshold <- rsq_threshold
  
  #### parametric fitting constraints
  # in the order: a0, ainf, ac50, hill
  constr.lower <- constr_lower
  constr.upper <- constr_upper

  x <- df$conc_um
  y <- df$readout_value_pctAC
  # add a tiny bit of noise to avoid problems with data without variance
  y <- y + rnorm(length(y), mean=0, sd=0.01)
  

  
  ### data preparation for fitting ...
  
  # sort by x ...
  sorted <- sort(x, index.return=TRUE)
  
  x  <- x[sorted$ix]
  y  <- y[sorted$ix]	
  
  logx <- log10(x) 
  
  n.conc <- length(unique(x))
  n.data <- length(x)
  n.rep  <- n.data / n.conc      # average # of replicates at each conc

  c.min <- min(x)
  c.max <- max(x)

  
  ### constant fit
  c.res <- constant.fit.results(logx, y, method="robust")

  ### np fit for parameter approximation
  np.res <- nonparametric.fit.results(logx, y, abs.acx.level=50)
  
  
  # for approximate curve shape diagnostics (sigmoid, more like bell-shaped ...) we need these quantities
  id.max     <- which.max(values(np.res$fit, logx))
  logx.max   <- logx[id.max]
  amax.data  <- values(np.res$fit, logx[id.max])

    # highest of the (possibly multiple) indexes where max fitted model value is reached
  max.id.max <- max(which(values(np.res$fit, logx)==amax.data))
  # amin,amax, r.sq are need for assessment of suitability of constant fit, further below.
  # If amax amin difference is small then we might revert to constant fitting, see further below.
  # The values can also be used for check for constant fit when finding bell-shaped curves.

  
  # set regression starting parameters from initial non-parametric fit results
  init.params <- c(np.res$a0, np.res$ainf, np.res$ac50, np.res$hill)
  
  # parametric fit, OLS
  # core function is already using try - if failed: $message = "fit failed" 
  p.res <- param.4PL.fit.results(logx, y, params=init.params, method="nls", abs.acx.level=abs.acx.level, constr.lower=constr.lower, constr.upper=constr.upper)
  
  
  ## ensure all 3 have same columns
  # np.res with column alpha, p.res with column rweights
  p.res$rweights <- NULL
  np.res$alpha <- NULL
  
  c.res$fit_method <- c.res$fit.type
  p.res$fit_method <- p.res$fit.type
  np.res$fit_method <- np.res$fit.type

  c.res$fit_type <- "constant"
  p.res$fit_type <- "param"
  np.res$fit_type <- "nonparam"
  
  c.res$fit.type <- NULL
  p.res$fit.type <- NULL
  np.res$fit.type <- NULL
  
  
  # calculate AC50 prioritization
  if (np.res$r.squared < rsq.threshold) {
    # bad fit of np -> constant curve
    c.res$fit_sel <- 1
    p.res$fit_sel <- 0
    np.res$fit_sel <- 0
  } else
    if ((np.res$ainf < (amax.data - 3*np.res$rse)) && (n.conc - (max.id.max/n.rep) >= 2)) {
      # bell shaped curve -> np fit
      c.res$fit_sel <- 0
      p.res$fit_sel <- 0
      np.res$fit_sel <- 1
    } else
      if ((p.res$r.squared < rsq.threshold) || (abs(p.res$amin - p.res$amax) < a.diff.thresh) || p.res$message %in% c("fit failed","initial par violates constraints")) {
        # bad fit of p OR amin - amax too low  -> constant fit
        c.res$fit_sel <- 1
        p.res$fit_sel <- 0
        np.res$fit_sel <- 0
      } else {
        # sigmoidal fit
        c.res$fit_sel <- 0
        p.res$fit_sel <- 1
        np.res$fit_sel <- 0
      } 
  
  

  all.res <- as.data.frame(rbind(c.res,np.res,p.res))
  
  if (best_fit) {all.res <- all.res %>% filter(fit_sel==1)}
  
  return(all.res %>% dplyr::select(-fit) %>% rowwise() %>% to_df(.))
}


# create data frame from results
to_df <- function(dat){
  
  dat$rse <- as.numeric(dat$rse)                   
  dat$r.squared <- as.numeric(dat$r.squared)          
  dat$a0  <- as.numeric(dat$a0)                   
  dat$ainf    <- as.numeric(dat$ainf)               
  dat$ac50   <- as.numeric(dat$ac50)                
  dat$hill   <- as.numeric(dat$hill)                
  dat$abs.acx   <- as.numeric(dat$abs.acx )             
  dat$abs.acx.level     <- as.numeric(dat$abs.acx.level)     
  dat$amax    <- as.numeric(dat$amax  )              
  dat$c.amax   <- as.numeric(dat$c.amax )              
  dat$amin     <- as.numeric(dat$amin )              
  dat$c.amin     <- as.numeric(dat$c.amin)            
  dat$convergence   <- as.character(dat$convergence  )         
  dat$message         <- as.character(dat$message  )      
  dat$fit_method    <- as.character(dat$fit_method) 
  dat$fit_type    <- as.character(dat$fit_type) 
  dat$fit_sel    <- as.character(dat$fit_sel) 
  # dat$fit   <- as.character(dat$fit) # results in many linebreaks in output file - dont do it
  dat$fit   <- NULL
  
  return(dat)
}

