
#' Compute the PCA
#' 
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
#' utils::data(Exp1_R25_prot, package='DAPARdata')
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), 'missing', level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op='>=', th=1)
#' obj <- MetaCellFiltering(obj, indices, cmd='delete')
#' res.pca <- wrapper.pca(obj$new)
#' 
#' 
#' @export
#' 
wrapper.pca <- function(obj, var.scaling=TRUE, ncp=NULL){
  
  if (! requireNamespace("FactoMineR", quietly = TRUE)) {
    stop("Please install FactoMineR: BiocManager::install('FactoMineR')")
  }
  
  # require(FactoMineR)
  if (missing(obj))
    stop("'obj' is required")
  else if (nrow(obj)==0)
    return(NULL)

  if (is.null(var.scaling)) 
    var.scaling <- TRUE
  
  res.pca <- NULL
  if (length(which(is.na(Biobase::exprs(obj)))) == 0){
    
  if (is.null(ncp)){
    nmax <- 12
    y <- Biobase::exprs(obj)
    nprot <- dim(y)[1]
    n <- dim(y)[2] # If too big, take the number of conditions.
    
    if (n > nmax){
      n <- length(unique(Biobase::pData(obj)$Condition))
    }
    
    
    ncp <- min(n, nmax)
  }
  # parameters available to the user
  variance.scaling <- TRUE
  
  res.pca <- FactoMineR::PCA(Biobase::exprs(obj), scale.unit = var.scaling, ncp=ncp, graph=FALSE)
  # si warning pour les missing values, le reproduire dans l'interface graphique
  }
  return(res.pca)
}



#' Plots the variables of PCA
#' 
#' 
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
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' res.pca <- wrapper.pca(Exp1_R25_pept)
#' plotPCA_Var(res.pca)
#' 
#' 
#' @export
#' 
plotPCA_Var <- function(res.pca, chosen.axes=c(1,2)){
  
  if (! requireNamespace("factoextra", quietly = TRUE)) {
    stop("Please install factoextra: BiocManager::install('factoextra')")
  }
  #plot.PCA(res.pca, choix="var", axes = chosen.axes, title="Sample factor map (PCA)")
  #require(factoextra)
  # Colorer en fonction du cos2: qualit? de repr?sentation
  if (is.null(res.pca)){return(NULL)}
  factoextra::fviz_pca_var(res.pca, axes = chosen.axes, col.var = "cos2",
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
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' res.pca <- wrapper.pca(Exp1_R25_pept)
#' plotPCA_Ind(res.pca)
#' 
#' 
#' @export
#' 
plotPCA_Ind <- function(res.pca, chosen.axes=c(1,2)){
  
  if (! requireNamespace("factoextra", quietly = TRUE)) {
    stop("Please install factoextra: BiocManager::install('factoextra')")
  }
  
  
  #plot.PCA(res.pca, choix="ind", axes = chosen.axes, select = 0:-1, title="Protein factor map (PCA)")
  if (is.null(res.pca))
    return(NULL)
  
  #require(factoextra)
  plot <- factoextra::fviz_pca_ind(res.pca,  axes = chosen.axes, geom="point")
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
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' res.pca <- wrapper.pca(Exp1_R25_pept, ncp=6)
#' plotPCA_Eigen(res.pca)
#' 
#' @export
#' 
plotPCA_Eigen <- function(res.pca){
  if (is.null(res.pca)){return(NULL)}
  eig.val <- res.pca$eig
  barplot(eig.val[, 2], 
          names.arg = 1:nrow(eig.val), 
          main = "Variances Explained by PCs (%)",
          xlab = "Principal Components",
          ylab = "Percentage of variances",
          col ="steelblue")
  # Add connected line segments to the plot
  lines(x = 1:nrow(eig.val), eig.val[, 2], 
        type = "b", pch = 19, col = "red")
  
}



#' @title Plots the eigen values of PCA with the highcharts library
#' 
#' @param res.pca xxx
#' 
#' @return A histogram
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' res.pca <- wrapper.pca(Exp1_R25_pept, ncp=6)
#' plotPCA_Eigen_hc(res.pca)
#' 
#' @import highcharter
#' 
#' @export
#' 
plotPCA_Eigen_hc <- function(res.pca){
  if (is.null(res.pca)){return(NULL)}
  hc <- highchart() %>%
    hc_yAxis_multiples(list(title = list(text = "% of variances"),lineWidth = 0, labels = list(format = "{value}%"), max = 100), 
                       list(title = list(text = "Cumulative % of variances"), opposite = FALSE, max = 100),
                       list(title = list(text = "Eigen values"), opposite = TRUE, labels = list(format = "{value}%") )) %>%
    hc_xAxis(title = "Principal Components", categories = rownames(res.pca$eig)) %>%
    hc_add_series(data.frame(y=res.pca$eig[,2]), type="column",  name = "% of variances", yAxis = 0) %>%
    hc_add_series(data.frame(y=res.pca$eig[,3]), type="line",  color="darkblue",name = "Cumulative % of variances", marker = "diamond", color = "#FF7900", yAxis = 0) %>%
    #hc_add_series(data.frame(y=res.pca$eig[,1]),  type="line",  lineWidth = 4,name = "Eigen values", color = "orange", yAxis = 2) %>%
    #hc_tooltip(crosshairs = TRUE, headerFormat = "<b>{point.x}</b><br>") %>%
    hc_legend(enabled = TRUE)
  # hc_plotOptions(column = list(colorByPoint = TRUE, colors = SiteOTD$Colors))
  
  
}
