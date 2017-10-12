

##' This function is a wrappper to the function groupGO from the
##' package clusterProfiler. Given a vector of genes/proteins, it returns the 
##' GO profile at a specific level. It returns a groupGOResult instance. 
##' 
##' 
##' @title Calculates the GO profile of a vector of genes/proteins at specific 
##' level
##' @param data A vector of ID (genes or proteins !!DIRE LESQUELS!!)
##' @param idFrom xxxxx
##' @param idTo xxxxx
##' @param orgdb xxxxx
##' @param ont xxxxx
##' @param level Specific GO Level
##' @param readable xxxxx
##' @return xxxxx
##' @author xxxxx
group_GO <- function(data, idFrom, idTo, orgdb, ont, level, readable=FALSE){
    
    require(as.character(orgdb),character.only = TRUE)
    gene <- bitr(data, fromType=idFrom, toType=idTo, OrgDb=orgdb)
    gene.id = gene$ENTREZID
    ggo <- groupGO(gene = gene.id, 
                 OrgDb = orgdb, 
                 ont = ont, 
                 level = level, 
                 readable= readable)
    
    return(ggo)
}


##' This function is a wrappper to the function enrichGO from the
##' package \code{\link{clusterProfiler}}. Given a vector of genes/proteins, it returns the 
##' GO xxxxxxx.  
##' 
##' @title Calculates the GO xxxxxxxxx
##' @param data A vector of ID (genes or proteins !!DIRE LESQUELS!!)
##' @param idFrom xxxxx
##' @param idTo xxxxx
##' @param orgdb xxxxx
##' @param ont One of "MF", "BP", and "CC" subontologies
##' @param readable xxxxx
##' @param pval The qvalue cutoff (same parameter as in the function \code{enrichGO} of the package \code{\link{clusterProfiler}})
##' @param universe xxxxxxxxx
##' @return A groupGOResult instance.
##' @author xxxxx
enrich_GO <- function(data, idFrom, idTo, orgdb, ont, readable=FALSE, pval, universe)
{
  tmp <- which(is.na(data))
  if (length(tmp) > 0){
      data <- data[-which(is.na(data))]
  }
  
    gene <- bitr(data, fromType=idFrom, toType=idTo, OrgDb=orgdb)
    if (is.null(gene)){return (NULL)}
    gene.id = gene$ENTREZID
    
    ego <- enrichGO(gene = gene.id, OrgDb = orgdb, ont = ont, 
                  pAdjustMethod="BH", 
                  pvalueCutoff=pval, 
                  readable=readable,
                  universe = NULL)   
    
    return(ego)
}



##' Returns the universe = totality of ID of an OrgDb annotation package
##' 
##' @title Returns the totality of ENTREZ ID (gene id) of an OrgDb annotation 
##' package. Careful : org.Pf.plasmo.db : no ENTREZID but ORF
##' @param orgdb a Bioconductor OrgDb annotation package 
##' @return A vector of ENTREZ ID (totality if the ID for the package) 
##' @author Florence Combes
univ_AnnotDbPkg <- function(orgdb){
    
    require(as.character(orgdb),character.only = TRUE)
    univ <- AnnotationDbi::keys(get(orgdb), keytype="ENTREZID")
    #different syntax for 'org.Pf.plasmo.db' package
    #univ<-keys(get(orgdb), keytype="ORF")
    return(univ)
}


##' This method returns a \code{\link{MSnSet}} object with the results
##' of the Gene Ontology analysis.
##' 
##' @title Returns a \code{\link{MSnSet}} object with the results of
##' the GO analysis performed with the \code{\link{clusterProfiler}} package. 
##' @param obj An object of the class \code{\link{MSnSet}}
##' @param ggo_res The object returned by the function \code{group_GO} of the package \code{DAPAR}
##' or the function \code{groupGO} of the package \code{\link{clusterProfiler}}
##' @param ego_res The object returned by the function \code{enrich_GO} of the package \code{DAPAR}
##' or the function \code{enrichGO} of the package \code{\link{clusterProfiler}}
##' @param organism The parameter xxx of the function
##' @param ontology One of "MF", "BP", and "CC" subontologies
##' @param levels A vector of the different GO grouping levels to save
##' @param pvalueCutoff The qvalue cutoff (same parameter as in the function \code{enrichGO} of the package \code{\link{clusterProfiler}})
##' @param typeUniverse  The type of background to be used. Values are 'Entire Organism', 'Entire dataset' or 'Custom'. In the latter
##' case, a file should be uploaded by the user
##' @return An object of the class \code{\link{MSnSet}}
##' @author Samuel Wieczorek
GOAnalysisSave <- function (obj, ggo_res=NULL, ego_res=NULL, organism, ontology, levels, pvalueCutoff, typeUniverse){
    if (is.null(ggo_res) && is.null(ego_res)){
        warning("Neither ggo or ego analysis has  been completed.")
        return(NULL)}
    
    if (!is.null(ggo_res)){
        text <- paste("Group analysis on ", organism)
        obj@processingData@processing <- c(obj@processingData@processing, text)

        obj@experimentData@other$GGO_analysis <- list(ggo_res = ggo_res,
                                                    organism = organism,
                                                    ontology = ontology,
                                                    levels = levels)
    }
    
    if (!is.null(ego_res)){
        text <- paste("Enrichment analysis on", organism)
        obj@processingData@processing <- c(obj@processingData@processing, text)
        
        obj@experimentData@other$EGO_analysis <- list(ego_res = ego_res,
                                                      organism = organism,
                                                      ontology = ontology,
                                                      PAdjustMethod = "BH",
                                                      pvalueCutoff = pvalueCutoff,
                                                      typeUniverse = typeUniverse)
    }
    
    
    return(obj)
}



##' A barplot of GO classification analysis
##' 
##' @title A barplot which shows the result of a GO classification, using the package \code{\link{highcharter}}
##' @param ggo The result of the GO classification, provides either by the function
##' \code{group_GO} in the package \code{DAPAR} or the function \code{groupGO} 
##' in the package \code{\link{clusterProfiler}}
##' @param maxRes An integer which is the maximum number of classes to display in the plot 
##' @param title The title of the plot
##' @return A barplot 
##' @author Samuel Wieczorek
barplotGroupGO_HC <- function(ggo, maxRes=5, title=""){
    
    dat <- ggo@result
    nRes <- min(maxRes, nrow(dat))
    
    n <- which(dat[,"Count"]==0)
    if (length(n) > 0){dat <- dat[-which(dat[,"Count"]==0),]}
    dat <- dat[order(dat[,"Count"], decreasing=TRUE),]
    dat <- dat[seq(1:nRes),]
    
        
    h1 <-  highchart() %>%
    my_hc_chart(chartType = "bar") %>%
    hc_title(text =title) %>%
    hc_add_series(dat[,"Count"]) %>%
    hc_legend(enabled = FALSE) %>%
    #hc_colors(myColors) %>%
        hc_tooltip(enabled = FALSE) %>%
     hc_xAxis(categories = dat[,"Description"], title = list(text = ""))

return(h1)
}



##' A barplot of GO enrichment analysis
##' 
##' @title A barplot that shows the result of a GO enrichment, using the package \code{\link{highcharter}}
##' @param ego The result of the GO enrichment, provides either by the function
##' \code{enrichGO} in the package \code{DAPAR} or the function \code{enrichGO} 
##' of the package \code{\link{clusterProfiler}}
##' @param maxRes The maximum number of categories to display in the plot 
##' @param title The title of the plot
##' @return A barplot 
##' @author Samuel Wieczorek
barplotEnrichGO_HC <- function(ego, maxRes = 5, title=NULL){
   if (is.null(ego)){return(NULL)}
    dat <- ego@result
    nRes <- min(maxRes, nrow(dat))
    
    n <- which(dat[,"Count"]==0)
    if (length(n) > 0){dat <- dat[-which(dat[,"Count"]==0),]}
    dat <- dat[order(dat[,"pvalue"], decreasing=FALSE),]
    dat <- dat[seq(1:nRes),]
    
    
    colfunc <- colorRampPalette(c("red","royalblue"))
    nbBreaks <- 20*nRes
    pal <- colfunc(nbBreaks)
    t <- log(dat[,"pvalue"])
    d <- (max(t) - min(t))/nbBreaks
    base <- seq(from=min(t), to=max(t), by = d)
    myColorsIndex <- unlist(lapply(t, function(x){last(which(x > base))}))
    myColorsIndex[which(is.na(myColorsIndex))] <- 1
    myColors <- pal[myColorsIndex]
    
    dat[,"pvalue"] <- format(dat[,"pvalue"], digits=2)
    
    h1 <- highchart() %>%  
        hc_title(title = title) %>%
        hc_yAxis(title = list(text = "Count")) %>% 
        hc_xAxis(categories = dat[,"Description"]) %>% 
        hc_add_series(data = data.frame(y=dat[,"Count"]), type = "bar", 
                      dataLabels = list(enabled = FALSE),
                      colorByPoint = TRUE) %>%
        hc_colors(myColors) %>%
        hc_tooltip(headerFormat= '', 
                   pointFormat = "{point.y}") %>%
        my_hc_ExportMenu(filename = "GOEnrich_barplot") %>%
        hc_legend(enabled = FALSE) %>%
        hc_plotOptions(bar = list(
            pointWidth=60,
            dataLabels = list(enabled = TRUE)))
    
     return(h1)
}


##' A scatter plot of GO enrichment analysis
##' 
##' @title A dotplot that shows the result of a GO enrichment, using the package \code{\link{highcharter}}
##' @param ego The result of the GO enrichment, provides either by the function
##' enrichGO in \code{DAPAR} or the function \code{enrichGO} of the packaage \code{\link{clusterProfiler}}
##' @param maxRes The maximum number of categories to display in the plot
##' @param title The title of the plot
##' @return A dotplot 
##' @author Samuel Wieczorek
scatterplotEnrichGO_HC <- function(ego, maxRes = 10, title=NULL){
    
    dat <- ego@result
    nRes <- min(maxRes, nrow(dat))
    dat$GeneRatio <- unlist(lapply(dat$GeneRatio, function(x){
        as.numeric(unlist(strsplit(x,'/'))[1]) / as.numeric(unlist(strsplit(x,'/'))[2])
    }))
    
    
    n <- which(dat$GeneRatio==0)
    if (length(n) > 0){dat <- dat[-which(dat$GeneRatio==0),]}
    dat <- dat[order(dat$GeneRatio, decreasing=TRUE),]
    dat <- dat[seq(1:nRes),]
    
    
    colfunc <- colorRampPalette(c("red","royalblue"))
    nbColors <- 5
    
    pal <- colfunc(nbColors)
    t <- log(dat$p.adjust)
    d <- (max(t) - min(t))/nbColors
    base <- seq(from=min(t), to=max(t), by = d)
    tmpList <- lapply(t, function(x){
                                  if (x == min(t)){ ind <- 1}
                                  else {ind <- which(x > base)[length(which(x > base))]}
                                  
    })
    
    myColorsIndex <- unlist(tmpList)
    
    df <- data.frame(x=c(0:(nRes-1)),
                        y=dat$GeneRatio,
                        z=dat$Count,
                        color=pal[myColorsIndex],
                        colorSegment=pal[myColorsIndex],
                        pAdjust = format(dat$p.adjust, digits=2),
                        name = dat[,"Description"])
    
    
    txt_tooltip <- paste("<b> p.adjust </b>: {point.pAdjust} <br> ", 
                         "<b> Count </b>: {point.z} <br> ",
                         sep="")
    
    h1 <-  highchart() %>%
        hc_title(title = title) %>%
        my_hc_chart(chartType = "bubble") %>%
        hc_add_series(df) %>%
        hc_legend(enabled = FALSE) %>%
        hc_xAxis(type = "category", categories = df$name)  %>%
        hc_yAxis(title = list(text = "Gene Ratio"))  %>%
        hc_tooltip(headerFormat= '',
                   pointFormat = txt_tooltip) %>%
        my_hc_ExportMenu(filename = "GOEnrich_dotplot")
    

    
    return(h1)
}