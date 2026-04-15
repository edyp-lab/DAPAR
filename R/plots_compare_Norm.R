
#' @title Builds a plot from a dataframe
#' 
#' @description Wrapper to the function that plot to compare the quantitative 
#' proteomics data before and after normalization.
#'
#' @param objBefore A dataframe that contains quantitative data before
#' normalization.
#'
#' @param objAfter A dataframe that contains quantitative data after
#' normalization.
#'
#' @param condsForLegend A vector of the conditions (one condition
#' per sample).
#'
#' @param ... arguments for palette
#'
#' @return A plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' 
#' data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept
#' conds <- Biobase::pData(obj)[, "Condition"]
#' objAfter <- wrapper.normalizeD(
#' obj = obj, method = "QuantileCentering",
#' conds = conds, type = "within conditions"
#' )
#' wrapper.compareNormalizationD_HC(obj, objAfter, conds,
#' pal = ExtendPalette(2))
#'
#' @export
#'
#'
wrapper.compareNormalizationD_HC <- function(objBefore,
    objAfter,
    condsForLegend = NULL,
    ...) {
    qDataBefore <- Biobase::exprs(objBefore)
    qDataAfter <- Biobase::exprs(objAfter)

    compareNormalizationD_HC(qDataBefore, 
        qDataAfter, 
        conds = condsForLegend, 
        ...)
}




#' @title Builds a plot from a dataframe. Same as compareNormalizationD but
#' uses the library \code{plotly}
#' 
#' @description 
#' Plot to compare the quantitative proteomics data before and after
#' normalization using the package \code{plotly}
#'
#'
#' @param qDataBefore A dataframe that contains quantitative data before
#' normalization.
#'
#' @param qDataAfter A dataframe that contains quantitative data after
#' normalization.
#'
#' @param keyId xxx
#'
#' @param conds A vector of the conditions (one condition
#' per sample).
#'
#' @param pal xxx
#'
#' @param subset.view xxx
#'
#' @param n An integer that is equal to the maximum number of displayed points.
#' This number must be less or equal to the size
#' of the dataset. If it is less than it, it is a random selection
#'
#' @param type scatter or line
#'
#' @return A plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot
#' qDataBefore <- Biobase::exprs(obj)
#' conds <- Biobase::pData(obj)[, "Condition"]
#' id <- Biobase::fData(obj)[, 'Protein_IDs']
#' pal <- ExtendPalette(2)
#' objAfter <- wrapper.normalizeD(obj,
#' method = "QuantileCentering",
#' conds = conds, type = "within conditions"
#' )
#'
#' n <- 1
#' compareNormalizationD_HC(
#' qDataBefore = qDataBefore,
#' qDataAfter = Biobase::exprs(objAfter), 
#' keyId = id, 
#' pal = pal, 
#' n = n,
#' subset.view = seq_len(n),
#' conds = conds)
#'
#' @import plotly
#' @importFrom utils str
#'
#' @export
#'
compareNormalizationD_HC <- function(qDataBefore,
    qDataAfter,
    keyId = NULL,
    conds = NULL,
    pal = NULL,
    subset.view = NULL,
    n = 1,
    type = "scatter") {
    
    pkgs.require('RColorBrewer')
    
  if (is.null(conds)) {
    warning("'conds' is null.")
    return(NULL)
  }
  if (n < 0 || n > nrow(qDataBefore)){
    warning("'n' must be a positive integer not null and less than the total number
      of entities. Set to default value: 0.2")
    n <- ceiling(0.2 * nrow(qDataBefore))
  }
  
  if (is.null(keyId)) {
    keyId <- seq_len(length(qDataBefore))
  }
  
  if (!is.null(subset.view) && length(subset.view) > 0) {
    keyId <- keyId[subset.view]
    if (nrow(qDataBefore) > 1) {
      if (length(subset.view) == 1) {
        qDataBefore <- t(qDataBefore[subset.view, ])
        qDataAfter <- t(qDataAfter[subset.view, ])
      } else {
        qDataBefore <- qDataBefore[subset.view, ]
        qDataAfter <- qDataAfter[subset.view, ]
      }
      n <- 100
    }
  } else {
    subset.view <- seq_len(n)
  }
  
  if (!match(type, c("scatter", "line"))) {
    warning("'type' must be equal to 'scatter' or 'line'.")
    return(NULL)
  }
  
  # Truncate dataset
  ind <- sample(seq_len(nrow(qDataBefore)), min(n, length(subset.view)))
  keyId <- keyId[ind]
  if (nrow(qDataBefore) > 1) {
    if (length(ind) == 1) {
      qDataBefore <- t(qDataBefore[ind, ])
      qDataAfter <- t(qDataAfter[ind, ])
    } else {
      qDataBefore <- qDataBefore[ind, ]
      qDataAfter <- qDataAfter[ind, ]
    }
  }
  
  myColors <- NULL
  if (is.null(pal)) {
    warning("Color palette set to default.")
    myColors <- GetColorsForConditions(conds, 
                                       ExtendPalette(length(unique(conds))))
  } else if (length(pal) != length(unique(conds))) {
    warning("The color palette has not the same dimension as 
              the number of samples")
    myColors <- GetColorsForConditions(conds, 
                                       ExtendPalette(length(unique(conds))))
  } else {
    myColors <- GetColorsForConditions(conds, pal)
  }
  
  x <- qDataBefore
  y <- qDataAfter / qDataBefore
  
  ## Colors definition
  legendColor <- unique(myColors)
  txtLegend <- unique(conds)
  
  shapes <- c(
    "circle", "square", "diamond",
    "triangle-up", "triangle-down",
    "cross", "x"
  )
  
  p <- plot_ly()
  
  for (i in seq_along(conds)) {
    p <- p |> add_markers(
      x = x[, i],
      y = y[, i],
      name = colnames(x)[i],
      text = paste0("Id: ", keyId),
      hoverinfo = "text",
      marker = list(
        color = myColors[i],
        symbol = shapes[(i - 1) %% length(shapes) + 1],
        size = 8
      )
    )
  }
  
  p <- p |> layout(
    legend = list(
      orientation = "h", 
      x = 0, 
      y = -0.15, 
      xanchor = "left",
      yanchor = "top"
    )
  )
  
  p
}
