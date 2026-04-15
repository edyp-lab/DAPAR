


#' @title Builds a densityplot from a dataframe
#' 
#' @description 
#' Densityplot of quantitative proteomics data over samples.
#'
#'
#' @param obj xxx
#'
#' @param legend A vector of the conditions (one condition
#' per sample).
#'
#' @param pal xxx
#'
#' @return A density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' densityPlotD_HC(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' pal <- ExtendPalette(2, "Dark2")
#' densityPlotD_HC(Exp1_R25_pept, pal = pal)
#'
#' @import plotly
#'
#' @export
#'
densityPlotD_HC <- function(obj,
                            legend = NULL,
                            pal = NULL) {
    
    pkgs.require('stats')
    
    
    if (is.null(obj)) {
        warning("The dataset is NULL and cannot be shown")
        return(NULL)
    } else if (nrow(obj) == 0) {
        warning("The dataset is empty and cannot be shown")
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    conds <- Biobase::pData(obj)$Condition

    if (is.null(legend)) {
        legend <- Biobase::pData(obj)[, "Condition"]
    }

    myColors <- NULL
    if (is.null(pal)) {
        myColors <- GetColorsForConditions(conds, 
            ExtendPalette(length(unique(conds))))
    } else {
        if (length(pal) != length(unique(conds))) {
            warning("The color palette has not the same dimension as the 
                number of samples. Set to default.")
            myColors <- GetColorsForConditions(conds, 
                ExtendPalette(length(unique(conds))))
        } else {
            myColors <- GetColorsForConditions(conds, pal)
        }
    }

    p <- plotly::plot_ly()
    
    for (i in seq_len(ncol(qData))) {
      
      dens <- stats::density(qData[, i], na.rm = TRUE)
      
      p <- p |>
        plotly::add_trace(
          x = dens$x,
          y = dens$y,
          type = "scatter",
          mode = "lines",
          name = legend[i],
          line = list(color = myColors[i]),
          hovertemplate = paste0(
            "<b>", legend[i], "</b>: %{y:.2f}<extra></extra>"
          )
        )
    }
    
    p <- p |>
      plotly::layout(
        title = "Density plot",
        xaxis = list(title = "log(Intensity)"),
        yaxis = list(title = "Density"),
        margin = list(t = 60, b = 60),
        legend = list(
          orientation = "h",
          x = 0,
          y = -0.15,
          xanchor = "left",
          yanchor = "top"
        )
      )
    
    return(p)
}
