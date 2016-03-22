#' Test dataset
#'
#' @name test
#' @docType data
#' @keywords data
#' @description Partial (small) dataset for unit tests containing
#' missing values.
#' @format An object of class \code{\link{MSnSet}}
NULL


#' Test dataset
#'
#' @name testWithoutNA
#' @docType data
#' @keywords data
#' @description Partial (small) dataset for unit tests without any
#' missing values.
#' @format An object of class \code{\link{MSnSet}}
NULL


#' UPSpep25 dataset
#'
#' @name UPSpep25
#' @docType data
#' @keywords data
#' @description This dataset is the final outcome of a quantitative mass spectrometry-based proteomic 
#' analysis of two samples containing different concentrations of 48 human proteins (UPS1 standard 
#' from Sigma-Aldrich) within a constant yeast background (see Giai Gianetto et al. (2016) for details). 
#' It contains the abundance values of the different human and yeast peptides identified and quantified 
#' in these two conditions. The two conditions represent the measured abundances of peptides when 
#' respectively 25fmol and 10fmol of UPS1 human proteins were mixed with the yeast extract before mass 
#' spectrometry analyses. Three technical replicates were acquired for each condition.
#' 
#' To identify and quantify peptides, spectra were searched using MaxQuant (version 1.5.1.2) 
#' against the Uniprot database, the UPS database and the frequently observed contaminants database. 
#' Maximum false discovery rates were set to 0.01 by employing a reverse 
#' database strategy.
#' 
#'
#' The dataset is either available as a CSV file (see inst/extdata/UPSpep25.txt), or as a \code{\link{MSnSet}}
#' structure (UPSpep25). In the latter case, the quantitative data are those of the raw intensities.
#' @usage data(UPSpep25)
#' @format An object of class \code{\link{MSnSet}} related to peptide
#' quantification. It contains 6 samples divided into two conditions
#' (25fmol and 10fmol) and 13918 peptides.
#' 
#' The data frame exprs(UPSpep25) contains six columns that are the quantitation of peptides for the six replicates. 
#' 
#' The data frame fData(UPSpep25) contains the meta data about the peptides.
#' 
#' The data frame pData(UPSpep25) contains the experimental design and gives few informations about the samples.
#' 
#' @references Cox J., Hein M.Y., Luber C.A., Paron I., Nagaraj N., Mann M. 
#' Accurate proteome-wide label-free quantification by delayed normalization and 
#' maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep, 13(9):2513-26.
#' 
#' Giai Gianetto, Q., Combes, F., Ramus, C., Bruley, C., Coute, Y., Burger, T. (2016). 
#' Calibration plot for proteomics: A graphical tool to visually check the assumptions underlying 
#' FDR control in quantitative experiments. Proteomics, 16(1), 29-32.
#' 
#' @keywords datasets
NULL
