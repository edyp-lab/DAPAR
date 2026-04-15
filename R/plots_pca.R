

#' @title Compute the PCA
#'
#' @param obj xxx
#'
#' @param var.scaling The dimensions to plot
#'
#' @param ncp xxxx
#'
#' @return A xxxxxx
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
#' res.pca <- wrapper.pca(obj$new)
#'
#' @export
#'
wrapper.pca <- function(obj, var.scaling = TRUE, ncp = NULL) {
    
    pkgs.require('FactoMineR')

    # require(FactoMineR)
    if (missing(obj)) {
        stop("'obj' is required")
    } else if (nrow(obj) == 0) {
        return(NULL)
    }

    if (is.null(var.scaling)) {
        var.scaling <- TRUE
    }

    res.pca <- NULL
    if (length(which(is.na(Biobase::exprs(obj)))) == 0) {
        if (is.null(ncp)) {
            nmax <- 12
            y <- Biobase::exprs(obj)
            nprot <- dim(y)[1]
            n <- dim(y)[2] # If too big, take the number of conditions.

            if (n > nmax) {
                n <- length(unique(Biobase::pData(obj)$Condition))
            }


            ncp <- min(n, nmax)
        }
        # parameters available to the user
        variance.scaling <- TRUE

        res.pca <- FactoMineR::PCA(
            Biobase::exprs(obj), 
            scale.unit = var.scaling, 
            ncp = ncp, 
            graph = FALSE)
        # si warning pour les missing values, le reproduire 
        # dans l'interface graphique
    }
    return(res.pca)
}




#' @title Plots variables of PCA
#'
#' @param res.pca xxx
#'
#' @param chosen.axes The dimensions to plot
#'
#' @return A plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' res.pca <- wrapper.pca(Exp1_R25_pept)
#' plotPCA_Var(res.pca)
#'
#' @export
#'
plotPCA_Var <- function(res.pca, chosen.axes = c(1, 2)) {
    pkgs.require('factoextra')
    
    # plot.PCA(res.pca, choix="var", axes = chosen.axes, 
    # title="Sample factor map (PCA)")
    # require(factoextra)
    # Colorer en fonction du cos2: qualit? de repr?sentation
    if (is.null(res.pca)) {
        return(NULL)
    }
    factoextra::fviz_pca_var(
        res.pca,
        axes = chosen.axes, 
        col.var = "cos2",
        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
        repel = TRUE # ?vite le chevauchement de texte
    )
}



#' @title Plots individuals of PCA
#'
#' @param res.pca xxx
#'
#' @param chosen.axes The dimensions to plot
#'
#' @return A plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' res.pca <- wrapper.pca(Exp1_R25_pept)
#' plotPCA_Ind(res.pca)
#'
#' @export
#'
plotPCA_Ind <- function(res.pca, chosen.axes = c(1, 2)) {
    pkgs.require('factoextra')

    if (is.null(res.pca)) {
        return(NULL)
    }

    plot <- factoextra::fviz_pca_ind(
        res.pca, 
        axes = chosen.axes, 
        geom = "point"
        )
    plot
}



#' @title Plots the eigen values of PCA
#'
#' @param res.pca xxx
#'
#' @return A histogram
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' res.pca <- wrapper.pca(Exp1_R25_pept, ncp = 6)
#' plotPCA_Eigen(res.pca)
#'
#' @export
#'
plotPCA_Eigen <- function(res.pca) {
    
    pkgs.require('graphics')
    
    if (is.null(res.pca)) {
        return(NULL)
    }
    eig.val <- res.pca$eig
    graphics::barplot(
        eig.val[, 2],
        names.arg = seq_len(nrow(eig.val)),
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col = "steelblue"
    )
    # Add connected line segments to the plot
    graphics::lines(
        x = seq_len(nrow(eig.val)), 
        eig.val[, 2],
        type = "b", 
        pch = 19, 
        col = "red"
    )
}



#' @title Plots the eigen values of PCA with the plotly library
#'
#' @param res.pca xxx
#'
#' @return A histogram
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package='DAPARdata')
#' res.pca <- wrapper.pca(Exp1_R25_pept, ncp = 6)
#' plotPCA_Eigen_hc(res.pca)
#'
#' @import plotly
#'
#' @export
#'
plotPCA_Eigen_hc <- function(res.pca) {
    if (is.null(res.pca)) {
        return(NULL)
    }
  eig <- res.pca$eig
  pcs <- rownames(eig)
  
  var_pct <- eig[, 2]
  cum_var <- eig[, 3]
  eigenvalues <- eig[, 1]
  
  plotly::plot_ly() |>
    plotly::add_trace(
      x = pcs,
      y = var_pct,
      type = "bar",
      name = "% of variances",
      yaxis = "y1"
    ) |>
    plotly::add_trace(
      x = pcs,
      y = cum_var,
      type = "scatter",
      mode = "lines+markers",
      name = "Cumulative % of variances",
      line = list(color = "#FF7900"),
      marker = list(symbol = "diamond"),
      yaxis = "y1"
    ) |>
    plotly::layout(
      xaxis = list(title = "Principal Components"),
      yaxis = list(
        title = "% of variances",
        range = c(0, 100),
        ticksuffix = "%"),
      yaxis2 = list(
        title = "Eigen values",
        overlaying = "y",
        side = "right"),
      legend = list(orientation = "h")
    )
}
