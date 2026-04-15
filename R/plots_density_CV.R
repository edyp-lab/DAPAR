
#' @title Distribution of CV of entities
#' 
#' @description Builds a densityplot of the CV of entities in the 
#' Biobase::exprs() table. of an object \code{MSnSet}. The variance is 
#' calculated for each condition present in the dataset (see the slot 
#' \code{'Condition'} in the \code{Biobase::pData()} table).
#'
#' @param obj An object of class \code{MSnSet}
#'
#' @param ... arguments for palette.
#'
#' @return A density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' wrapper.CVDistD_HC(Exp1_R25_pept)
#'
#' @export
#'
wrapper.CVDistD_HC <- function(obj, ...) {
    if (nrow(obj) == 0) {
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    conds <- Biobase::pData(obj)[, "Condition"]
    CVDistD_HC(qData, conds, ...)
}



#'
#' @title Distribution of CV of entities
#' 
#' @description 
#' Builds a densityplot of the CV of entities in the Biobase::exprs() table
#' of a object. The CV is calculated for each condition present
#' in the dataset (see the slot \code{'Condition'} in the 
#' \code{Biobase::pData()} table)
#'
#' @param qData A dataframe that contains quantitative data.
#'
#' @param conds A vector of the conditions (one condition per sample).
#'
#' @param pal xxx
#'
#' @return A density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' conds <- Biobase::pData(Exp1_R25_pept)[, "Condition"]
#' CVDistD_HC(Biobase::exprs(Exp1_R25_pept), conds)
#' pal <- ExtendPalette(2, "Dark2")
#' CVDistD_HC(Biobase::exprs(Exp1_R25_pept), conds, pal)
#'
#' @import plotly
#'
#' @export
#'
CVDistD_HC <- function(qData,
    conds = NULL,
    pal = NULL) {
    
    pkgs.require('stats')
    
    
    if (is.null(conds)) {
        warning("The vector of conditions is empty. The plot cannot be drawn.")
        return(NULL)
    }

    conditions <- unique(conds)
    n <- length(conditions)


    if (is.null(pal)) {
        pal <- ExtendPalette(n)
    } else {
        if (length(pal) != n) {
            warning("The color palette has not the same dimension as the 
                number of samples. Set to default.")
            pal <- ExtendPalette(n)
        }
    }
    
    p <- plotly::plot_ly()

    minX <- Inf
    maxX <- -Inf
    
    for (i in seq_len(n)) {
      
      idx <- which(conds == conditions[i])
      
      if (length(idx) > 1) {
        t <- apply(
          qData[, idx, drop = FALSE], 1,
          function(x) {
            m <- mean(x, na.rm = TRUE)
            if (is.na(m) || m == 0) return(NA)
            100 * stats::var(x, na.rm = TRUE) / m
          }
        )
        
        t <- t[!is.na(t)]
        
        if (length(t) > 1) {
          dens <- stats::density(t)
          
          minX <- min(minX, dens$x)
          xmaxY <- dens$x[which.max(dens$y)]
          maxX <- max(maxX, 10 * (xmaxY - minX))
          
          p <- p |>
            plotly::add_trace(
              x = dens$x,
              y = dens$y,
              type = "scatter",
              mode = "lines",
              name = conditions[i],
              line = list(color = pal[i]),
              hovertemplate = paste0(
                "<b>", conditions[i], "</b>: %{y:.2f}<extra></extra>"
              )
            )
        }
      }
    }
    
    if (!is.finite(minX) || !is.finite(maxX)) {
      minX <- NULL
      maxX <- NULL
    }
    
    p <- p |>
      plotly::layout(
        xaxis = list(
          title = "CV(log(Intensity))",
          range = if (!is.null(minX)) c(minX, maxX) else NULL, 
          zeroline = FALSE
        ),
        yaxis = list(title = "Density"),
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
