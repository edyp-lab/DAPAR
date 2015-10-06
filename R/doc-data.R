#' UPSprotx2 dataset
#'
#' @name UPSprotx2
#' @docType data
#' @keywords data
#' @description This dataset results from a controlled relative quantification
#' proteomics experiment where the commercial Sigma mix UPS1 human proteins
#' were spiked in a similar yeast lysate in 2 different concentrations
#' (with a ratio of 2).
#' As a consequence, it can be used to benchmark the quality of a
#' statistical analysis: in the ideal case, after the differential analysis,
#' only and all the human proteins should have thus been selected.
#' The dataset is either available as a CSV file
#' (see inst/extdata/proteinGroups-UPSx2.txt), or as a \code{\link{MSnSet}}
#' structure (\code{data(UPSprotx2)}). In the latter case, the quantitative
#' data are those of the raw intensities.
#' @format An object of class \code{\link{MSnSet}} related to proteins
#' quantification. It contains 6 samples divided into two conditions
#' (5fmol and 10fmol) and 2394 proteins.
NULL

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
