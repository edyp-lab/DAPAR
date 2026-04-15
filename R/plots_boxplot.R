
#' @title Builds a boxplot from a dataframe using the package \code{plotly}
#' @param obj Numeric matrix
#' @param conds xxx
#' @param keyId xxxx
#' @param legend A vector of the conditions (one condition per sample).
#' @param pal A basis palette for the boxes which length must be equal
#' to the number of unique conditions in the dataset.
#' @param subset.view A vector of index indicating which rows to highlight
#' @return A boxplot
#' @author Samuel Wieczorek, Anais Courtier, Enora Fremy
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot
#' conds <- legend <- Biobase::pData(obj)$Condition
#' key <- "Protein_IDs"
#' pal <- ExtendPalette(length(unique(conds)))
#' boxPlotD_HC(obj, conds, key, legend, pal, seq_len(10))
#' @import plotly
#' @export
boxPlotD_HC <- function(
    obj,
    conds,
    keyId = NULL,
    legend = NULL,
    pal = NULL,
    subset.view = NULL) {
    pkgs.require(c('stats', "grDevices", "RColorBrewer"))
    

    if (is.null(obj)) {
        warning("The dataset is NULL and cannot be shown")
        return(NULL)
    } else if (nrow(obj) == 0) {
        warning("The dataset is empty and cannot be shown")
        return(NULL)
    } else {
        qData <- Biobase::exprs(obj)
        samples <- colnames(qData)
    }


    pkgs.require('plotly')

    if (missing(conds)) {
        stop("'conds' is missing.")
    }

    if (length(subset.view) == 0) {
        subset.view <- NULL
    }

    if (is.null(legend)) {
        legend <- conds
        for (i in unique(conds)) {
            legend[which(conds == i)] <- paste0(i, "_",
                seq_len(length(which(conds == i))))
        }
    }
    
    conds <- as.factor(conds)
    cond_levels <- levels(conds)

    myColors <- NULL
    if (is.null(pal)) {
      warning("Color palette set to default.")
      myColors <- GetColorsForConditions(unique(conds),
                                         ExtendPalette(length(unique(conds))))
    } else {
      if (length(pal) != length(unique(conds))) {
        warning("The color palette has not the same dimension as
                the number of conditions")
        myColors <- GetColorsForConditions(unique(conds),
                                           ExtendPalette(length(unique(conds))))
      } else {
        myColors <- GetColorsForConditions(unique(conds), pal)
      }
    }
    
    cond_color_map <- stats::setNames(myColors, unique(conds))
    
    p <- plotly::plot_ly()
    
    for (cl in cond_levels) {
      
      idx <- which(conds == cl)
      
      df_cond <- data.frame(
        sample = rep(samples[idx], each = nrow(qData)),
        value = as.vector(qData[, idx])
      )
      
      p <- p |> plotly::add_trace(
        data = df_cond,
        x = ~sample,
        y = ~value,
        type = "box",
        name = cl,
        fillcolor = cond_color_map[[cl]],
        marker = list(color = cond_color_map[[cl]]),
        line = list(color = "black"),
        boxpoints = "outliers",
        showlegend = FALSE
      )
    }
    
    if (!is.null(subset.view)) {
      
      if (is.null(keyId)) stop("'keyId' is missing.")
      if (!keyId %in% colnames(Biobase::fData(obj))) {
        stop("'keyId' does not belong to metadata")
      }
      
      pal2 <- ExtendPalette(length(subset.view), "Dark2")
      
      for (i in seq_along(subset.view)) {
        
        idx <- subset.view[i]
        
        overlay_df <- data.frame(
          sample = samples,
          value = as.numeric(qData[idx, ])
        )
        
        p <- p |> plotly::add_trace(
          data = overlay_df,
          x = ~sample,
          y = ~value,
          type = "scatter",
          mode = "lines+markers",
          line = list(color = pal2[i], dash = "dot"),
          marker = list(color = pal2[i], size = 6),
          showlegend = FALSE
        )
      }
    }
    
    p <- p |> plotly::layout(
      xaxis = list(title = "Samples"),
      yaxis = list(title = "Log (intensity)"),
      margin = list(b = 60),
      showlegend = FALSE
    )
    
    return(p)
}
