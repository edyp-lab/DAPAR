

#' This function is a wrappper to the function groupGO from the
#' package `clusterProfiler`. Given a vector of genes/proteins,
#' it returns the GO profile at a specific level. It returns a groupGOResult
#' instance.
#'
#' @title Calculates the GO profile of a vector of genes/proteins at a given
#' level of the Gene Ontology
#'
#' @param data A vector of ID (among ENSEMBL, ENTREZID, GENENAME, REFSEQ,
#' UNIGENE, UNIPROT -can be different according to organisms)
#'
#' @param idFrom character indicating the input ID format (among ENSEMBL,
#' ENTREZID, GENENAME, REFSEQ, UNIGENE, UNIPROT)
#'
#' @param orgdb annotation Bioconductor package to use (character format)
#'
#' @param ont on which ontology to perform the analysis (MF, BP or CC)
#'
#' @param level level of the ontolofy to perform the analysis
#'
#' @param readable TRUE or FALSE (default FALSE)
#'
#' @return GO profile at a specific level
#'
#' @author Florence Combes
#'
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(10)]
#' if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
#' stop("Please install org.Sc.sgd.db: 
#'             BiocManager::install('org.Sc.sgd.db')")
#' }
#' library(org.Sc.sgd.db)
#' ggo <- group_GO(
#'     data = Biobase::fData(obj)$Protein.IDs, idFrom = "UNIPROT",
#'     orgdb = "org.Sc.sgd.db", ont = "MF", level = 2
#' )
#'
#' @export
#'
group_GO <- function(data, idFrom, orgdb, ont, level, readable = FALSE) {
    
    pkg.orgdb <- as.character(orgdb)

    pkgs.require(c('clusterProfiler', as.character(orgdb)))
    
    if (idFrom == "UNIPROT") {
        gene <- clusterProfiler::bitr(
            data, 
            fromType = idFrom, 
            toType = "ENTREZID", 
            OrgDb = orgdb)
        if (is.null(gene)) {
            return(NULL)
        }
        gene.id <- gene$ENTREZID
    } else {
        gene.id <- data
    }

    ggo <- clusterProfiler::groupGO(
        gene = gene.id,
        OrgDb = orgdb,
        ont = ont,
        level = level,
        readable = readable
    )

    return(ggo)
}


#' @title Calculates GO enrichment classes for a given list of proteins/genes 
#' ID. It results an enrichResult instance.
#' 
#' @description 
#' This function is a wrappper to the function enrichGO from the package 
#' `clusterProfiler`. Given a vector of genes/proteins, it returns an 
#' enrichResult instance.
#'
#' @param data A vector of ID (among ENSEMBL, ENTREZID, GENENAME, REFSEQ,
#' UNIGENE, UNIPROT -can be different according to organisms)
#'
#' @param idFrom character indicating the input ID format (among ENSEMBL,
#' ENTREZID, GENENAME, REFSEQ, UNIGENE, UNIPROT)
#'
#' @param orgdb annotation Bioconductor package to use (character format)
#'
#' @param ont One of "MF", "BP", and "CC" subontologies
#'
#' @param readable TRUE or FALSE (default FALSE)
#'
#' @param pval The qvalue cutoff (same parameter as in the function
#' \code{enrichGO} of the package `clusterProfiler`)
#'
#' @param universe a list of ID to be considered as the background for
#' enrichment calculation
#'
#' @return A groupGOResult instance.
#'
#' @author Florence Combes
#'
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(10)]
#' if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
#' stop("Please install org.Sc.sgd.db: 
#'             BiocManager::install('org.Sc.sgd.db')")
#' }
#' library(org.Sc.sgd.db)
#' univ <- univ_AnnotDbPkg("org.Sc.sgd.db") # univ is the background

#' ego <- enrich_GO(
#'     data = Biobase::fData(obj)$Protein.IDs, idFrom = "UNIPROT",
#'     orgdb = "org.Sc.sgd.db", ont = "MF", pval = 0.05, universe = univ
#' )
#'
#' @export
#'
enrich_GO <- function(
    data, 
    idFrom, 
    orgdb, 
    ont, 
    readable = FALSE, 
    pval, 
    universe) {
    
    pkgs.require('clusterProfiler')

    tmp <- which(is.na(data))
    if (length(tmp) > 0) {
        data <- data[-which(is.na(data))]
    }


    if (idFrom == "UNIPROT") {
        gene <- clusterProfiler::bitr(
            data, 
            fromType = idFrom, 
            toType = "ENTREZID", 
            OrgDb = orgdb)
        if (is.null(gene)) {
            return(NULL)
        }
        gene.id <- gene$ENTREZID
    } else {
        gene.id <- data
    }

    ego <- clusterProfiler::enrichGO(
        gene = gene.id, OrgDb = orgdb, ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = pval,
        readable = readable,
        universe = NULL
    )

    return(ego)
}




#'
#' @title Returns the totality of ENTREZ ID (gene id) of an OrgDb annotation
#' package. Careful : org.Pf.plasmo.db : no ENTREZID but ORF
#'
#' @description
#' Function to compute the `universe` argument for the \code{enrich_GO}
#' function, in case this latter should be the entire organism.
#' Returns all the ID of the OrgDb annotation package for the corresponding
#' organism.
#' 
#' @param orgdb a Bioconductor OrgDb annotation package
#'
#' @return A vector of ENTREZ ID
#'
#' @author Florence Combes
#'
#' @export
#' 
#' @examples
#' if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
#' stop("Please install org.Sc.sgd.db: 
#'             BiocManager::install('org.Sc.sgd.db')")
#' }
#' library(org.Sc.sgd.db)
#' univ_AnnotDbPkg("org.Sc.sgd.db")
#'
univ_AnnotDbPkg <- function(orgdb) {
    
    pkgs.require(c('AnnotationDbi', as.character(orgdb)))
    
    
    univ <- AnnotationDbi::keys(get(orgdb), keytype = "ENTREZID")
    # different syntax for 'org.Pf.plasmo.db' package
    # univ<-keys(get(orgdb), keytype="ORF")
    return(univ)
}




#' @title Returns an \code{MSnSet} object with the results of the GO analysis 
#' performed with the functions \code{enrichGO} and/or \code{groupGO} of the 
#' `clusterProfiler` package.
#' 
#' @description 
#' This method returns an \code{MSnSet} object with the results of the Gene 
#' Ontology analysis.
#'
#' @param obj An object of the class \code{MSnSet}
#'
#' @param ggo_res The object returned by the function \code{group_GO} of the
#' package \code{DAPAR} or the function \code{groupGO} of the package
#'  `clusterProfiler`
#'
#' @param ego_res The object returned by the function \code{enrich_GO} of the
#' package \code{DAPAR} or the function \code{enrichGO} of the package
#' `clusterProfiler`
#'
#' @param organism The parameter OrgDb of the functions \code{bitr},
#' \code{groupGO} and \code{enrichGO}
#'
#' @param ontology One of "MF", "BP", and "CC" subontologies
#'
#' @param levels A vector of the different GO grouping levels to save
#'
#' @param pvalueCutoff The qvalue cutoff (same parameter as in the function
#' \code{enrichGO} of the package `clusterProfiler`)
#'
#' @param typeUniverse  The type of background to be used. Values are
#' 'Entire Organism', 'Entire dataset' or 'Custom'. In the latter case, a file
#' should be uploaded by the user
#'
#' @return An object of the class \code{MSnSet}
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @examples 
#' NULL
#'
GOAnalysisSave <- function(obj,
    ggo_res = NULL,
    ego_res = NULL,
    organism,
    ontology,
    levels,
    pvalueCutoff,
    typeUniverse) {
    if (is.null(ggo_res) && is.null(ego_res)) {
        warning("Neither ggo or ego analysis has  been completed.")
        return(NULL)
    }

    if (!is.null(ggo_res)) {
        text <- paste("Group analysis on ", organism)
        
        obj@experimentData@other$GGO_analysis <- list(
            ggo_res = ggo_res,
            organism = organism,
            ontology = ontology,
            levels = levels
        )
    }

    if (!is.null(ego_res)) {
        text <- paste("Enrichment analysis on", organism)
        
        obj@experimentData@other$EGO_analysis <- list(
            ego_res = ego_res,
            organism = organism,
            ontology = ontology,
            PAdjustMethod = "BH",
            pvalueCutoff = pvalueCutoff,
            typeUniverse = typeUniverse
        )
    }


    return(obj)
}



#' @title A barplot which shows the result of a GO classification, using the
#' package \code{plotly}
#'
#' @param ggo The result of the GO classification, provides either by the
#' function \code{group_GO} in the package \code{DAPAR} or the function
#' \code{groupGO} in the package `clusterProfiler`
#'
#' @param maxRes An integer which is the maximum number of classes to display
#' in the plot
#'
#' @param title The title of the plot
#'
#' @return A barplot
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(10)]
#' if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
#' stop("Please install org.Sc.sgd.db: 
#'             BiocManager::install('org.Sc.sgd.db')")
#' }
#' library(org.Sc.sgd.db)
#' univ <- univ_AnnotDbPkg("org.Sc.sgd.db")
#' ggo <- group_GO(
#'     data = Biobase::fData(obj)$Protein.IDs, idFrom = "UNIPROT",
#'     orgdb = "org.Sc.sgd.db", ont = "MF", level = 2
#' )
#' barplotGroupGO_HC(ggo)
#'
barplotGroupGO_HC <- function(ggo, maxRes = 5, title = "") {
    dat <- ggo@result
    nRes <- min(maxRes, nrow(dat))

    dat <- dat[dat[, "Count"] != 0, ]
    dat <- dat[order(dat[, "Count"], decreasing = TRUE), ]
    dat <- dat[seq_len(min(nRes, nrow(dat))), ]

    dat$Description <- factor(dat$Description, levels = rev(dat$Description))
    
    p <- plotly::plot_ly(
      data = dat,
      x = ~Count,
      y = ~Description,
      type = "bar",
      orientation = "h"
    ) |>
      plotly::layout(
        title = title,
        xaxis = list(title = ""),
        yaxis = list(title = ""),
        showlegend = FALSE
      )
    
    return(p)
}


#' A barplot of GO enrichment analysis
#'
#' @title A barplot that shows the result of a GO enrichment, using the
#' package \code{plotly}
#'
#' @param ego The result of the GO enrichment, provides either by the function
#' \code{enrichGO} in the package \code{DAPAR} or the function \code{enrichGO}
#' of the package `clusterProfiler`
#'
#' @param maxRes The maximum number of categories to display in the plot
#'
#' @param title The title of the plot
#'
#' @return A barplot
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(10)]
#' if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
#' stop("Please install org.Sc.sgd.db: 
#'             BiocManager::install('org.Sc.sgd.db')")
#' }
#' library(org.Sc.sgd.db)
#' univ <- univ_AnnotDbPkg("org.Sc.sgd.db")
#' ego <- enrich_GO(
#'     data = Biobase::fData(obj)$Protein.IDs, idFrom = "UNIPROT",
#'     orgdb = "org.Sc.sgd.db", ont = "MF", pval = 0.05, universe = univ
#' )
#' barplotEnrichGO_HC(ego)
#'
barplotEnrichGO_HC <- function(ego, maxRes = 5, title = NULL) {
    
    pkgs.require('dplyr')
    
    if (is.null(ego)) {
        return(NULL)
    }
    dat <- ego@result
    nRes <- min(maxRes, nrow(dat))

    dat <- dat[dat$Count != 0, ]
    dat <- dat[order(dat[, "pvalue"], decreasing = FALSE), ]
    dat <- dat[seq_len(min(nRes, nrow(dat))), ]


    colfunc <- grDevices::colorRampPalette(c("red", "royalblue"))
    nbBreaks <- 20 * nRes
    pal <- colfunc(nbBreaks)
    t <- log(dat[, "pvalue"])
    d <- (max(t) - min(t)) / nbBreaks
    base <- seq(from = min(t), to = max(t), by = d)
    myColorsIndex <- unlist(lapply(t, function(x) {
        dplyr::last(which(x > base))
    }))
    myColorsIndex[is.na(myColorsIndex)] <- 1
    myColors <- pal[myColorsIndex]

    hover_text <- paste0(
      "<b>Description:</b> ", dat$Description, "<br>",
      "<b>Count:</b> ", dat$Count, "<br>",
      "<b>pvalue:</b> ", format(dat$pvalue, digits = 2)
    )
    
    p <- plotly::plot_ly(
      data = dat,
      x = ~reorder(Description, Count),
      y = ~Count,
      type = "bar",
      marker = list(color = myColors),
      text = hover_text,
      hoverinfo = "text"
    ) |>
      layout(
        title = title,
        xaxis = list(title = ""),
        yaxis = list(title = "Count")
      ) |>
      plotly::config(displayModeBar = TRUE)
    
    return(p)
}


#' A scatter plot of GO enrichment analysis
#'
#' @title A dotplot that shows the result of a GO enrichment, using the
#' package \code{plotly}
#'
#' @param ego The result of the GO enrichment, provides either by the function
#' enrichGO in \code{DAPAR} or the function \code{enrichGO} of the packaage
#' `clusterProfiler`
#'
#' @param maxRes The maximum number of categories to display in the plot
#'
#' @param title The title of the plot
#'
#' @return A dotplot
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot
#' if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
#' stop("Please install org.Sc.sgd.db: 
#'             BiocManager::install('org.Sc.sgd.db')")
#' }
#' library(org.Sc.sgd.db)
#' univ <- univ_AnnotDbPkg("org.Sc.sgd.db")
#' ego <- enrich_GO(
#'     data = Biobase::fData(obj)$Protein.IDs, idFrom = "UNIPROT",
#'     orgdb = "org.Sc.sgd.db", ont = "MF", pval = 0.05, universe = univ
#' )
#' scatterplotEnrichGO_HC(ego)

scatterplotEnrichGO_HC <- function(ego, maxRes = 10, title = NULL) {
  
  if (is.null(ego))
    return(NULL)
  
  pkgs.require('grDevices')
  
  dat <- ego@result
  nRes <- min(maxRes, nrow(dat))
  
  dat$GeneRatio <- vapply(dat$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }, numeric(1))
  
  dat <- dat[dat$GeneRatio != 0, ]
  
  dat <- dat[order(dat$GeneRatio, decreasing = TRUE), ]
  dat <- dat[seq_len(min(nRes, nrow(dat))), ]
  
  colfunc <- grDevices::colorRampPalette(c("red", "royalblue"))
  nbColors <- 5
  pal <- colfunc(nbColors)
  
  t <- log(dat$p.adjust)
  d <- (max(t) - min(t)) / nbColors
  base <- seq(min(t), max(t), by = d)
  
  myColorsIndex <- vapply(t, function(x) {
    if (x == min(t)) {
      1
    } else {
      which(x > base)[length(which(x > base))]
    }
  }, numeric(1))
  
  colors <- pal[myColorsIndex]
  
  df <- data.frame(
    name = dat$Description,
    GeneRatio = dat$GeneRatio,
    Count = dat$Count,
    pAdjust = dat$p.adjust,
    color = colors,
    stringsAsFactors = FALSE
  )
  
  p <- plotly::plot_ly(
    data = df,
    x = ~name,
    y = ~GeneRatio,
    type = "scatter",
    mode = "markers",
    marker = list(
      color = ~color,
      size = ~Count * 2,
      opacity = 0.8
    ),
    text = ~paste0(
      "<b>", name, "</b><br>",
      "p.adjust: ", format(pAdjust, digits = 2), "<br>",
      "Count: ", Count
    ),
    hoverinfo = "text"
  ) |>
    plotly::layout(
      title = title,
      xaxis = list(title = ""),
      yaxis = list(title = "Gene Ratio")
    )
  
  return(p)
}