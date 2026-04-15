#' @title Displays a correlation matrix of the quantitative data of the
#' \code{Biobase::exprs()} table
#' 
#' @description Builds a correlation matrix based on a \code{MSnSet} object.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param rate A float that defines the gradient of colors.
#'
#' @param showValues xxx
#'
#' @return A colored correlation matrix
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' wrapper.corrMatrixD_HC(Exp1_R25_pept)
#'
#'
#' @export
#'
wrapper.corrMatrixD_HC <- function(obj, rate = 0.5, showValues = TRUE) {
    
    pkgs.require('stats')
    
    
    if (is.null(obj)) {
        warning("The dataset is NULL and cannot be shown")
        return(NULL)
    } else if (nrow(obj) == 0) {
        warning("The dataset is empty and cannot be shown")
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)
    data <- stats::cor(qData, use = "pairwise.complete.obs")
    corrMatrixD_HC(data, samplesData, rate, showValues)
}




#' @title Displays a correlation matrix of the quantitative data of the
#' \code{Biobase::exprs()} table.
#'
#' @param object The result of the \code{cor} function.
#'
#' @param samplesData A dataframe in which lines correspond to samples and
#' columns to the meta-data for those samples.
#'
#' @param rate The rate parameter to control the exponential law for
#' the gradient of colors
#'
#' @param showValues xxx
#'
#' @return A colored correlation matrix
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' samplesData <- Biobase::pData(Exp1_R25_pept)
#' res <- cor(qData, use = "pairwise.complete.obs")
#' corrMatrixD_HC(res, samplesData)
#'
#' @import plotly
#'
#' @export
#'
corrMatrixD_HC <- function(object, 
    samplesData = NULL, 
    rate = 0.5, 
    showValues = TRUE) {
    
    pkgs.require(c('stats', "dplyr", "tidyr", "tibble"))
    
    df <- as.data.frame(object)
    .sData <- samplesData
    if (!is.null(.sData)) {
        for (j in seq_len(ncol(df))) {
            names(df)[j] <- paste(as.character(.sData[j, 2:ncol(.sData)]),
                collapse = " "
            )
        }
    }
    is.num <- vapply(df, is.numeric, FUN.VALUE = NA)
    df[is.num] <- lapply(df[is.num], round, 2)
    mat <- as.matrix(df)
    labels <- colnames(mat)
    
    text_mat <- if (showValues) {
      matrix(sprintf("%.2f", mat), nrow = nrow(mat))
    } else {
      NULL
    }
    
    plotly::plot_ly(
      x = labels,
      y = labels,
      z = mat,
      type = "heatmap",
      colorscale = list(
        list(0, "#FF5733"),
        list(0.5, "#F8F5F5"),
        list(1, "#2E86C1")
      ),
      zmin = rate,
      zmax = 1,
      text = text_mat,
      texttemplate = if (showValues) "%{text}" else NULL,
      hovertemplate = paste(
        "%{y} ~ %{x}: <b>%{z:.2f}</b><extra></extra>"
      )
    ) |>
      plotly::layout(
        xaxis = list(title = "", side = "top"),
        yaxis = list(title = ""),
        margin = list(l = 100, r = 100)
      )
}
