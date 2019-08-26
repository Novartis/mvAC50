
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
library(matrixStats)


#########################################
# Functions to quantify signature effects
#########################################


maha.dist <- function(df, ref){
  ref.medians <- matrixStats::colMedians(as.matrix(ref), na.rm=TRUE)
  ref.cov <- cov(ref, use="complete.obs")
  maha <- mahalanobis(df, ref.medians, ref.cov,tol=1e-40)
  return(maha)
}

euc.dist <- function(df, ref) {
  ref.medians <- matrixStats::colMedians(as.matrix(ref), na.rm=TRUE)
  euc <- apply(df, 1, function(x) sqrt(sum((x - ref.medians) ^ 2)))
  return(euc)
}

cor.pearson <- function(df, ref) {
  ref.medians <- matrixStats::colMedians(as.matrix(ref), na.rm=TRUE)
  cor.p <- apply(df, 1, function(x) cor(x, ref.medians, use="na.or.complete"))
  return(cor.p)
}

cor.spearman <- function(df, ref) {
  ref.medians <- matrixStats::colMedians(as.matrix(ref), na.rm=TRUE)
  cor.s <- apply(df, 1, function(x) cor(x, ref.medians, method="spearman", use="na.or.complete"))
  return(cor.s)
}

dot.product <- function(df, ref) {
  ref.medians <- matrixStats::colMedians(as.matrix(ref), na.rm=TRUE)
  dot.p <- apply(df, 1, function(x) sum(x*ref.medians))
  return(dot.p)
}

vec.norm <- function(df) {
  vec.n <- apply(df, 1, function(x) sqrt(sum(x*x)))
  return(vec.n)
}

cosine <- function(df, ref) {
  ref.medians <- matrixStats::colMedians(as.matrix(ref), na.rm=TRUE)
  cosi <- apply(df, 1, function(x) sum(x*ref.medians) / ( sqrt(sum(x * x)) * sqrt(sum(ref.medians * ref.medians)) ) )
  return(cosi)
}

scalar.projection.AC <- function(df, ref) {
  ref.medians <- matrixStats::colMedians(as.matrix(ref), na.rm=TRUE)
  nr <- apply(df, 1, function(x) sqrt(sum((x * ( sum(x*ref.medians) / ( sqrt(sum(x * x)) * sqrt(sum(ref.medians * ref.medians)) ) ))^2 )))
  return(nr)
}

significance.weight <- function(df, threshold=3) {
  significance.w <- apply(df, 1, function(x) min(1,mean(abs(x)) / abs(threshold)))
  return(significance.w)
}

cosine.weight <- function(df, ref, threshold=3) {
  ref.medians <- matrixStats::colMedians(as.matrix(ref), na.rm=TRUE)
  cosi <- apply(df, 1, function(x) min(1,mean(abs(x)) / abs(threshold)) * (sum(x*ref.medians) / ( sqrt(sum(x * x)) * sqrt(sum(ref.medians * ref.medians))) ) )
  return(cosi)
}

num.readouts.changed <- function(df, threshold=3){
  num.readouts <- apply(df, 1, function(x) sum(abs(x) > abs(threshold)))
  return(num.readouts)
}


##########################
# wrapper around functions
##########################


#' Signature scores
#' 
#' signature_scores calculates a set of summary methods for the effect of a treatment on multivariate readouts. 
#' Summary methode quantify effects based on the signature data alone, or relative to a neutral or active control signature.
#' 
#' 
#' @param df wide dataframe with one row per treatment and multivariate readouts in columns
#' @param readouts list of columns containing multivariate readouts to be summarize, e.g. c("col1", "col2", ...)
#' @param keep_cols list of columns with metadata to report in results, e.g. c("col1", "col2", ...)
#' @param well_type column specifying samples as "SA", neutral controls as "NC", and active controls as "AC"
#' @param mv_methods list of multivariate summary methods to compute. Options are:
#' \itemize{
#'  \item "maha_AC" mahalanobis distance to the AC,
#'  \item "maha_NC" mahalanobis distance to the NC,
#'  \item "euc_AC" euclidean distance to the AC,r
#'  \item "euc_NC" euclidean distance to the NC,r
#'  \item "cor_p_AC" Pearson correlation to the AC,
#'  \item "cor_s_AC" Spearman correlation to the AC,
#'  \item "dot_p_AC" dot product to the AC,
#'  \item "cos_AC" cosine to the AC,
#'  \item "vec_norm" vector norm,
#'  \item "scalar_projection_AC" scalar projection to the AC,
#'  \item "cos_weight_AC" cosine to the AC multiplied by max(1, mean(abs(readout_values)) / threshold)
#'  \item "num_readouts_changed" number of readouts with abs(values) above threshold.
#'  }
#'  Default are all methods except the mahalanobis distances.
#' @param threshold threshold used in cos_weight_AC and num_readouts_changed. Default is 3, assuming the readout values are normalized as zscores or rscores, so that a threshold of 3 means 3 standard deviations.
#' @param keep_n_most_responsive_readouts number of the original readouts to keep which are most responsive to the AC condition
#' @param keep_readouts list of columns of the original readouts to keep for results, e.g. c("col1", "col2", ...)
#'
#' @return Dataframe with keep_cols, potential original readouts, and multivariate readouts with
#' readout_value the raw value
#' readout_value_pctAC the readout value normalized as percentage of the AC
#' readout_value_pctACNC the readout value normalized to NC and as percentage of the AC  
#' @export
#'
#' @examples
#' see vignettes
#' 
signature_scores <- function(df, readouts, keep_cols, well_type, mv_methods=c(), threshold=3, keep_n_most_responsive_readouts=0, keep_readouts=c()){

  mv_methods_most <- c("euc_AC","euc_NC","cor_p_AC","cor_s_AC","dot_p_AC","cos_AC","vec_norm","scalar_projection_AC","cos_weight_AC","num_readouts_changed")
  # mv_methods_default <- c("cos_weight_AC","scalar_projection_AC", "vec_norm","num_readouts_changed")
  # mv_methods_all <- c("maha_AC","maha_NC","euc_AC","euc_NC","cor_p_AC","cor_s_AC","dot_p_AC","cos_AC","vec_norm","scalar_projection_AC","cos_weight_AC","num_readouts_changed")
  
  # if well_types col with a name != "well_type" add a col for easier selections
  if (well_type != "well_type"){df$well_type <- df[[well_type]]}

  if (length(mv_methods)==0){mv_methods <- mv_methods_most}
  if (!("AC" %in% df$well_type & "NC" %in% df$well_type  & "SA" %in% df$well_type)) {
      print("pleaese define well_types: AC for active control, NC for neutral control and SA for sample")
  }
  else{
    df.ac <- df[df$well_type=="AC",(colnames(df) %in% readouts)]
    df.nc <- df[df$well_type=="NC",(colnames(df) %in% readouts)]
    df.readouts <- df[,(colnames(df) %in% readouts)]

    # identify x most responding genes - and univariate readouts to keep
    readouts_most_resp <- colnames(df.ac)[rank(-abs(matrixStats::colMedians(as.matrix(df.ac), na.rm=TRUE)))<=keep_n_most_responsive_readouts]
    keep_readouts <- unique(c(keep_readouts, readouts_most_resp, "well_type"))
    
    # generate results 
    df.res <- df[,(colnames(df) %in% c(keep_cols, keep_readouts))]
    if("maha_AC" %in% mv_methods) {df.res$maha_AC <- maha.dist(df.readouts, df.ac)}
    if("maha_NC" %in% mv_methods) {df.res$maha_NC <- maha.dist(df.readouts, df.nc)}
    if("euc_AC" %in% mv_methods) {df.res$euc_AC <- euc.dist(df.readouts, df.ac)}
    if("euc_NC" %in% mv_methods) {df.res$euc_NC <- euc.dist(df.readouts, df.nc)}
    if("cor_p_AC" %in% mv_methods) {df.res$cor_p_AC <- cor.pearson(df.readouts, df.ac)}
    if("cor_s_AC" %in% mv_methods) {df.res$cor_s_AC <- cor.spearman(df.readouts, df.ac)}
    if("dot_p_AC" %in% mv_methods) {df.res$dot_p_AC <- dot.product(df.readouts, df.ac)}
    if("cos_AC" %in% mv_methods) {df.res$cos_AC <- cosine(df.readouts, df.ac)}
    if("vec_norm" %in% mv_methods) {df.res$vec_norm <- vec.norm(df.readouts)}
    if("scalar_projection_AC" %in% mv_methods) {
      # needs cos_AC for calculation - can we get that into the real function?
      if (! "cos_AC" %in% mv_methods) {df.res$cos_AC <- cosine(df.readouts, df.ac)}
      df.res$scalar_projection_AC <- sign(df.res$cos_AC) * scalar.projection.AC(df.readouts, df.ac) # important to multiply with sign of the cosine - otherwise always positiv
      if (! "cos_AC" %in% mv_methods) {df.res$cos_AC <- NULL}
    }
    # df.res$significance_weight <- significance.weight(df.readouts)}
    if("cos_weight_AC" %in% mv_methods) {df.res$cos_weight_AC <- cosine.weight(df.readouts, df.ac, threshold)}
    if("num_readouts_changed" %in% mv_methods) {df.res$num_readouts_changed <- num.readouts.changed(df.readouts, threshold)}
    
    # make long version of df
    df.res.l <- df.res %>% gather(readout_name, readout_value, -keep_cols)
    df.res.l <- df.res.l %>% mutate(readout_type=ifelse(readout_name %in% keep_readouts, "univariate", "multivariate"))
    df.res.l <- df.res.l %>% mutate(readout_class=ifelse(readout_name %in% keep_readouts, "univariate",
                                                  ifelse(readout_name %in% c("cor_p_AC","cor_s_AC","cos_AC","cos_weight_AC"), "direction",
                                                  ifelse(readout_name %in% c("dot_p_AC","scalar_projection_AC"),"directionMagnitude",
                                                  ifelse(readout_name %in% c("vec_norm","euc_NC","maha_NC","num_readouts_changed"),"magnitude","AC_similarity")))))

    
    ## calculate pct AC signal based on grouping
    df.res.l <- inner_join(df.res.l, df.res.l %>% filter(well_type=="AC") %>% group_by(readout_name) %>% summarize(ac.median=median(readout_value)), by="readout_name")
    df.res.l <- inner_join(df.res.l, df.res.l %>% filter(well_type=="NC") %>% group_by(readout_name) %>% summarize(nc.median=median(readout_value)), by="readout_name")
    df.res.l <- df.res.l  %>% mutate( readout_value_pctAC=readout_value / ac.median * 100,
                                      readout_value_pctACNC=(readout_value-nc.median) / (ac.median-nc.median) * 100) %>%
                                      dplyr::select(-ac.median, -nc.median)
    
    return(df.res.l)  
  }
}


