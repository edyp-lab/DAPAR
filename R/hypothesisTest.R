
#'
#' @title Density plots of logFC values
#' 
#' @description 
#' This function show the density plots of Fold Change (the same as calculated
#' by limma) for a list of the comparisons of conditions in a differential
#' analysis.
#'
#' @param df_logFC A dataframe that contains the logFC values
#'
#' @param threshold_LogFC The threshold on log(Fold Change) to
#' distinguish between differential and non-differential data
#'
#' @param pal xxx
#'
#' @return A plotly density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(100)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(GetMetacell(obj), c("Missing POV", "Missing MEC"), level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' qData <- Biobase::exprs(obj$new)
#' sTab <- Biobase::pData(obj$new)
#' res <- limmaCompleteTest(qData, sTab, comp.type = "OnevsAll")
#' pal <- ExtendPalette(2, "Dark2")
#' hc_logFC_DensityPlot(res$logFC, threshold_LogFC = 1, pal = pal)
#'
#' @export
#' 
#' @import plotly
#'
hc_logFC_DensityPlot <- function(df_logFC,
                                 threshold_LogFC = 0,
                                 pal = NULL) {
    
    pkgs.require(c("stats", "RColorBrewer", "grDevices"))
    
  if (threshold_LogFC < 0) {
    warning("The parameter 'threshold_LogFC' must be positive or equal to zero.")
    return(NULL)
  }
  
  if (is.null(df_logFC) || ncol(df_logFC) == 0) {
    return(NULL)
  }
  
  if (is.null(pal)) {
    warning("Color palette set to default.")
    pal <- ExtendPalette(ncol(df_logFC), "Paired")
  } else if (length(pal) != ncol(df_logFC)) {
    warning("The color palette has not the same dimension as the 
              number of samples")
    pal <- ExtendPalette(ncol(df_logFC), "Paired")
  }
  
  nValues <- nrow(df_logFC) * ncol(df_logFC)
  nInf <- sum(df_logFC <= -threshold_LogFC)
  nSup <- sum(df_logFC >=  threshold_LogFC)
  nInside <- sum(abs(df_logFC) < threshold_LogFC)
  
  
  p <- plot_ly()
  
  maxY.inf <- 0
  maxY.inside <- 0
  maxY.sup <- 0
  minX <- Inf
  maxX <- -Inf
  
  
  for (i in seq_len(ncol(df_logFC))) {
    
    tmp <- stats::density(df_logFC[, i], na.rm = TRUE)
    
    minX <- min(minX, tmp$x)
    maxX <- max(maxX, tmp$x)
    
    maxY.inf <- max(maxY.inf, max(tmp$y[tmp$x <= -threshold_LogFC], 0))
    maxY.inside <- max(maxY.inside, max(tmp$y[tmp$x > -threshold_LogFC & tmp$x < threshold_LogFC], 0))
    maxY.sup <- max(maxY.sup, max(tmp$y[tmp$x >=  threshold_LogFC], 0))
    
    p <- p |> add_lines(
      x = tmp$x,
      y = tmp$y,
      name = colnames(df_logFC)[i],
      line = list(color = pal[i]),
      hovertemplate = paste0("<b>", colnames(df_logFC)[i], "</b><br>",
                             "y: %{y:.2f}<extra></extra>"),
      showlegend = TRUE
    )
  }
  
  p <- p |> plotly::layout(
    title = "log(FC) repartition",
    margin = list(t = 60, b = 60),
    xaxis = list(title = "log(FC)"),
    yaxis = list(title = "Density"),
    legend = list(
      orientation = "h", 
      x = 0, 
      y = -0.15, 
      xanchor = "left",
      yanchor = "top"
    ),
    shapes = list(
      list(
        type = "rect",
        x0 = -threshold_LogFC,
        x1 = threshold_LogFC,
        y0 = 0,
        y1 = 1,
        xref = "x",
        yref = "paper",
        fillcolor = "lightgrey",
        opacity = 0.5,
        line = list(width = 0)
      )
    )
  )
  
  if (threshold_LogFC > 0) {
    p <- p |> add_annotations(
      x = 10,
      y = maxY.inside-0.1,
      text = sprintf("n Filtered out = %d<br>(%.2f%%)", nInside, 100*nInside/nValues),
      showarrow = FALSE,
      arrowhead = 2,
      ax = 40,
      ay = -40,
      font = list(size = 18)
    )
  }
  
  if (threshold_LogFC >= minX) {
    p <- p |> add_annotations(
      x = mean(c(minX, -threshold_LogFC)),
      y = maxY.inf+0.1,
      text = sprintf("nInf = %d<br>(%.2f%%)", nInf, 100*nInf/nValues),
      showarrow = FALSE,
      font = list(color = "blue",
                  size = 18)
    )
  }
  
  if (threshold_LogFC <= maxX) {
    p <- p |> add_annotations(
      x = mean(c(maxX, threshold_LogFC)),
      y = maxY.sup+0.1,
      text = sprintf("nSup = %d<br>(%.2f%%)", nSup, 100*nSup/nValues),
      showarrow = FALSE,
      font = list(color = "blue",
                  size = 18)
    )
  }
  
  return(p)
}
