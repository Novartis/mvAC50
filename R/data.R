
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


#' Dose-dependent changes in gene expression for THP1 cells treated with beta agonists
#'
#' A dataset containing the concentration dependent changes in expression 
#' of 14 genes and a housekeeper gene upon treatment with compounds at 
#' 8 different concentrations
#'
#' @format A data frame with 24192 rows and 13 variables:
#' \describe{
#'   \item{row}{row on plate}
#'   \item{col}{column on plate}
#'   \item{assay_plate}{assay plate}
#'   \item{well_type}{classify wells into one of "SA" for sample treatment, "AC" for active control, and "NC" for neutral control}
#'   \item{conc_um}{compound concentration in uM}
#'   \item{inchi_key}{inchi_key for compounds}
#'   \item{smiles}{smiles encoding of compound structure}
#'   \item{biol_replicate}{biological replicates}
#'   \item{gene}{gene measured}
#'   \item{rscore}{rscore: number of robust standard deviations of gene expression away from neutral controls}
#'   \item{rscore_HKnorm}{rscore normalized to housekeeing gene }
#'   \item{cpd_name}{name of compound}
#'   \item{readout_set}{the 14 genes were measured in two sets of 7 genes + housekeeper}
#' }
#' @source NIBR
"thp1_gene_rscores"


#' Dose-dependent changes in cAMP levels for THP1 cells treated with beta agonists
#'
#' A dataset containing the concentration dependent changes in cAMP levels
#' of THP1 cells treated with compounds.
#'
#' @format A data frame with 24192 rows and 13 variables:
#' \describe{
#'   \item{row}{row on plate}
#'   \item{col}{column on plate}
#'   \item{assay_plate}{assay plate}
#'   \item{well_type}{classify wells into one of "SA" for sample treatment, "AC" for active control, and "NC" for neutral control}
#'   \item{conc_um}{compound concentration in uM}
#'   \item{inchi_key}{inchi_key for compounds}
#'   \item{smiles}{smiles encoding of compound structure}
#'   \item{biol_replicate}{biological replicates}
#'   \item{readout_name}{name of readout}
#'   \item{readout_value}{readout value from assay}
#'   \item{readout_value_pctAC}{activity normalized to percent of the active control}
#'   \item{readout_type}{type of readout}
#'   \item{cpd_name}{name of compound}
#'   \item{readout_set}{annotation for comparison with gene signature results}
#' }
#' @source NIBR
"thp1_cAMP"

                
#' Dose-dependent changes in gene expression of an EGFR signature in MCF10A cells upon treatment with EGFR inhibitors
#'
#' A dataset containing the concentration dependent changes in expression 
#' of the genes of an EGFR signature upon treatment with EGFR inhibitors at 
#' 6 different concentrations. The data is a subset of the L1000 data set
#' More details on the full dataset in Subramanian et al., Cell, 2017 (10.1016/j.cell.2017.10.049)
#'
#' @format A data frame with 610 rows and 392 variables:
#' \describe{
#'   \item{well_type}{classify wells into one of "SA" for sample treatment, "AC" for active control, and "NC" for neutral control}
#'   \item{gene_signature}{annotation of gene signature for comparison of datasets}
#'   \item{cpd_id}{standardized compound name}
#'   \item{cl_id}{standardized cell line name}
#'   \item{id}{original id from L1000}
#'   \item{SM_Dose}{compound concentration in uM}
#'   \item{CL_Center_Specific_ID}{original cell lines id}
#'   \item{SM_Name}{name of compounds}
#'   \item{SM_Pert_Type}{pertubation type as annotated in L1000}
#'   \item{SM_Time}{incubation time}
#'   \item{det_plate}{assay plate}
#'   \item{det_well}{plate well id}
#'   \item{...}{many columns named after genes containing changes in normalized gene expression}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138}
"lincs_egfr_sig"


#' Dose-dependent changes in gene expression of an TPX2 signature in MCF10A cells upon treatment with EGFR inhibitors
#'
#' A dataset containing the concentration dependent changes in expression 
#' of the genes of an TPX2 signature upon treatment with EGFR inhibitors at 
#' 6 different concentrations. The data is a subset of the L1000 data set
#' More details on the full dataset in Subramanian et al., Cell, 2017 (10.1016/j.cell.2017.10.049)
#'
#' @format A data frame with 24192 rows and 13 variables:
#' \describe{
#'   \item{well_type}{classify wells into one of "SA" for sample treatment, "AC" for active control, and "NC" for neutral control}
#'   \item{gene_signature}{annotation of gene signature for comparison of datasets}
#'   \item{cpd_id}{standardized compound name}
#'   \item{cl_id}{standardized cell line name}
#'   \item{id}{original id from L1000}
#'   \item{SM_Dose}{compound concentration in uM}
#'   \item{CL_Center_Specific_ID}{original cell lines id}
#'   \item{SM_Name}{name of compounds}
#'   \item{SM_Pert_Type}{pertubation type as annotated in L1000}
#'   \item{SM_Time}{incubation time}
#'   \item{det_plate}{assay plate}
#'   \item{det_well}{plate well id}
#'   \item{...}{many columns named after genes containing changes in normalized gene expression}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138}
"lincs_tpx2_sig"


#' GR50 growth inhibition data for EGFR inhibitors in MCF10A cells
#'
#' A dataset containing GR50s activities of EGFR inhibitors in MCF10A cells 
#' The data and algorithm is described in more details in Hafner et al 2016 (doi: 10.1038/nmeth.3853)
#'
#' @format A data frame with 24192 rows and 13 variables:
#' \describe{
#'   \item{GR50.cellLine}{cell line for which data was measured}
#'   \item{GR50.smallMolecule}{compound}
#'   \item{GR50}{GR50 as defined in the Hafner paper}
#'   \item{GRmax}{GRmax as defined in the Hafner paper}
#'   \item{GR_AOC}{GR_AOC as defined in the Hafner paper}
#'   \item{GEC50}{GEC50 as defined in the Hafner paper}
#'   \item{GRinf}{GRinf as defined in the Hafner paper}
#'   \item{GR50.Hill}{GR50.Hill as defined in the Hafner paper}
#'   \item{GR50.r2}{GR50.r2 as defined in the Hafner paper}
#'   \item{GR50.pval}{GR50.pval as defined in the Hafner paper}
#'   \item{cpd_id}{standardized compound name}
#'   \item{cl_id}{standardized cell line name}
#'   \item{GR50.e}{estimated GR50s. For cpds with flat curves, more extreme GR50 potencies were estimated than just the maximal concentrations.}
#' }
#' @source \url{http://lincs.hms.harvard.edu/db/datasets/20252/results}
"lincs_GR50s"











