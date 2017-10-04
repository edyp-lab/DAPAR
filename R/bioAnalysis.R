

##' This function is a wrappper to the function groupGO from the
##' package clusterProfiler. Given a vector of genes/proteins, it returns the 
##' GO profile at a specific level. It returns a groupGOResult instance. 
##' 
##' 
##' @title Calculates the GO profile of a vector of genes/proteins at specific 
##' level
##' @param data A vector of ID (genes or proteins !!DIRE LESQUELS!!)
##' @param idFrom
##' @param idTo
##' @param threshold_LogFC The threshold on log(Fold Change) to
##' distinguish between differential and non-differential data 
##' @param pi0Method The parameter pi0.method of the method adjust.p 
##' in the package \code{cp4p}
##' @return The computed FDR value (floating number)
##' @author Florence Combes
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- wrapper.mvImputation(Exp1_R25_pept[1:1000], "QRILC")
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' qData <- Biobase::exprs(obj)
##' samplesData <- Biobase::pData(obj)
##' labels <- Biobase::pData(obj)[,"Label"]
##' limma <- diffAnaLimma(qData,samplesData, labels, condition1, condition2)
##' diffAnaComputeFDR(limma)
group_GO <- function(data, idFrom, idTo, orgdb, ont, level, readable=TRUE){
    
    require(as.character(orgdb),character.only = TRUE)
    gene <- bitr(data, fromType=idFrom, toType=idTo, OrgDb=orgdb)
    gene.id = gene$ENTREZID
    ggo <- groupGO(gene = gene.id, 
                 OrgDb = orgdb, 
                 ont = ont, 
                 level = level, 
                 readable= TRUE)
    
    return(ggo)
}


##' This function is a wrappper to the function groupGO from the
##' package clusterProfiler. Given a vector of genes/proteins, it returns the 
##' GO profile at a specific level. It returns a groupGOResult instance. 
##' 
##' @title Calculates the GO profile of a vector of genes/proteins at specific 
##' level
##' @param data A vector of ID (genes or proteins !!DIRE LESQUELS!!)
##' @param threshold_PVal The threshold on p-pvalue to
##' distinguish between differential and non-differential data 
##' @param threshold_LogFC The threshold on log(Fold Change) to
##' distinguish between differential and non-differential data 
##' @param pi0Method The parameter pi0.method of the method adjust.p 
##' in the package \code{cp4p}
##' @return The computed FDR value (floating number)
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- wrapper.mvImputation(Exp1_R25_pept[1:1000], "QRILC")
##' condition1 <- '25fmol'
##' condition2 <- '10fmol'
##' qData <- Biobase::exprs(obj)
##' samplesData <- Biobase::pData(obj)
##' labels <- Biobase::pData(obj)[,"Label"]
##' limma <- diffAnaLimma(qData,samplesData, labels, condition1, condition2)
##' diffAnaComputeFDR(limma)
enrich_GO <- function(data, idFrom, idTo, orgdb, ont, readable=TRUE, pAdj, pval, universe)
{
  data <- data[-which(is.na(data))]
  
    ## ENRICHMENT : GO over-representation test
    
    gene <- bitr(data, fromType=idFrom, toType=idTo, OrgDb=orgdb)
    if (is.null(gene)){return (NULL)}
    gene.id = gene$ENTREZID
    
    ego<-enrichGO(gene = gene.id, OrgDb = orgdb, ont = ont, 
                  pAdjustMethod=pAdj, 
                  pvalueCutoff=pval, 
                  readable=TRUE,
                  universe = NULL)   
    
    return(ego)
}


##' xxxxxxxxxxxxxx
##' 
##' @title xxxxxxxxxxxxxx
##' @param c xxxxxxxxxxxxxx 
##' @return xxxxxxxxxxxxxx 
##' @author Samuel Wieczorek
getUniprotID_FromString <- function(x){
    x <- unlist(x)
    uniprotSepIndices <- which(x=="|")
    if (length(uniprotSepIndices) >= 2) {
        x <- paste0(x, collapse="")
        res <- substr(x,uniprotSepIndices[1], uniprotSepIndices[2]-2)

    } else { res <- NA}
    return(res)
}

##' xxxxxxxxxxxxxx
##' 
##' @title xxxxxxxxxxxxxx
##' @param dat xxxxxxxxxxxxxx 
##' @return xxxxxxxxxxxxxx 
##' @author Samuel Wieczorek
getUniprotID_FromVector <- function(dat){
    require(stringr)
    d <- str_split(dat, "|", Inf)
    uniprotID <- lapply(d,getUniprotID_FromString)

    return(unlist(uniprotID))
}


##' Returns the universe = totality of ID of an OrgDb annotation package
##' 
##' @title Returns the totality of ENTREZ ID (gene id) of an OrgDb annotation 
##' package
##' @param orgdb a Bioconductor OrgDb annotation package 
##' @return A vector of ENTREZ ID (totality if the ID for the package) 
##' @author Florence Combes
##' @examples 
##' 
##' Careful : org.Pf.plasmo.db : no 'ENTREZID' but 'ORF'
univ_AnnotDbPkg <- function(orgdb){
    
    require(as.character(orgdb),character.only = TRUE)
    univ <- keys(get(orgdb), keytype="ENTREZID")
    #different syntax for 'org.Pf.plasmo.db' package
    #univ<-keys(get(orgdb), keytype="ORF")
    return(univ)
}


##' xxxxxxxxxxxxxx
##' 
##' @title xxxxxxxxxxxxxx
##' @param obj An object of the class \class{MSnset} 
##' @param ggo_res The object returned by the function \code{group_GO} of the package DAPAR
##' or the function \code{groupGO} of the package \CRANpkg{clusterProfiler}
##' @param ego_res The object returned by the function \code{enrich_GO} of the package DAPAR
##' or the function \code{enrichGO} of the package \CRANpkg{clusterProfiler}
##' @param organism The parameter xxx of the function
##' @param ontology xxxxxxxxxxxxxx
##' @param level xxxxxxxxxxxxxx
##' @param PAdjustMethod xxxxxxxxxxxxxx
##' @param pvalueCutoff xxxxxxxxx
##' @param typeUniverse 
##' @return An object of the class \xxx{MSnset}
##' @author Samuel Wieczorek
##' @examples
GOAnalysisSave <- function (obj, ggo_res, ego_res, organism, ontology, level, PAdjustMethod, pvalueCutoff, typeUniverse){
    if (is.null(ggo_res) && is.null(ego_res)){
        warning("Neither ggo or ego analysis has  been completed.")
        return(NULL)}
    
    if (!is.null(ggo_res)){
        text <- paste("Group analysis on ", organism)
        obj@processingData@processing <- c(obj@processingData@processing, text)

        obj@experimentData@other$GGO_analysis <- list(ggo_res = ggo_res,
                                                    organism = organism,
                                                    ontology = ontology,
                                                    level = level)
    }
    
    if (!is.null(ego_res)){
        text <- paste("Enrichment analysis on", organism)
        obj@processingData@processing <- c(obj@processingData@processing, text)
        
        obj@experimentData@other$EGO_analysis <- list(ego_res = ego_res,
                                                      organism = organism,
                                                      ontology = ontology,
                                                      PAdjustMethod = PAdjustMethod,
                                                      pvalueCutoff = pvalueCutoff,
                                                      typeUniverse = typeUniverse)
    }
    
    
    return(obj)
}



##' A barplot of GO classification analysis
##' 
##' @title A barplot that shows the result of a GO classification, using the package \CRANpkg{highcharter}
##' @param ego The result of the GO classification, provides either by the function
##' \code{group_GO} in the package \xxxpkg{DAPAR} or the function \code{xxx} of the packaage xxxx
##' @param maxRes An integer which is the maximum number of classes to display in the plot 
##' @return A barplot 
##' @author Samuel Wieczorek
##' @examples 
barplotGroupGO_HC <- function(ggo, maxRes=5){
    
    dat <- ggo@result
    nRes <- min(maxRes, nrow(dat))
    
    n <- which(dat$Count==0)
    if (length(n) > 0){dat <- dat[-which(dat$Count==0),]}
    dat <- dat[order(dat$Count, decreasing=TRUE),]
    dat <- dat[seq(1:nRes),]
    
        
    h1 <-  highchart() %>%
    hc_chart(type = "bar") %>%
    hc_add_series(dat$Count) %>%
    hc_legend(enabled = FALSE) %>%
    #hc_colors(myColors) %>%
     hc_xAxis(categories = dat$Description, title = list(text = ""))

return(h1)
}

##' A barplot of GO enrichment analysis
##' 
##' @title A barplot that shows the result of a GO enrichment, using the package \CRANpkg{highcharter}
##' @param ego The result of the GO enrichment, provides either by the function
##' \code {enrichGO} in the package DAPAR or the function \code{xxx} of the packaage xxxx
##' @param maxRes An integer which is the maximum number of categories to display in the plot 
##' @return A barplot 
##' @author Samuel Wieczorek
##' @examples

barplotEnrichGO_HC <- function(ego, maxRes = 5){
   
    dat <- ego@result
    nRes <- min(maxRes, nrow(dat))
    
    n <- which(dat$Count==0)
    if (length(n) > 0){dat <- dat[-which(dat$Count==0),]}
    dat <- dat[order(dat$pvalue, decreasing=FALSE),]
    dat <- dat[seq(1:nRes),]
    
    
    colfunc <- colorRampPalette(c("red","royalblue"))
    nbBreaks <- 20*nRes
    pal <- colfunc(nbBreaks)
    t <- log(dat$pvalue)
    d <- (max(t) - min(t))/nbBreaks
    base <- seq(from=min(t), to=max(t), by = d)
    myColorsIndex <- unlist(lapply(t, function(x){last(which(x > base))}))
    myColorsIndex[which(is.na(myColorsIndex))] <- 1
    myColors <- pal[myColorsIndex]
    
    dat$pvalue <- format(dat$pvalue, digits=2)
    
    h1 <- highchart() %>%  
        hc_yAxis(title = list(text = "Count")) %>% 
        hc_xAxis(categories = dat$Description) %>% 
        hc_add_series(data = dat, type = "bar", hcaes(x = Description, y = Count),
                      dataLabels = list(enabled = TRUE, format='pval={point.pvalue}'),
                      colorByPoint = TRUE) %>%
        hc_colors(myColors) %>%
        my_hc_ExportMenu(filename = "GOEnrich_barplot") %>%
        hc_legend(enabled = FALSE) %>%
        hc_plotOptions(bar = list(
            pointWidth=60,
            dataLabels = list(enabled = TRUE)))
    
     return(h1)
}


##' xxxxxxxxxxxxxx
##' 
##' @title A dotplot that shows the result of a GO enrichment, using the package \CRANpkg{highcharter}
##' @param ego The result of the GO enrichment, provides either by the function
##' enrichGO in DAPAR or bye the function \code{xxx} of the packaage xxxx
##' maxRes An integer which is the maximum number of categories to display in the plot 
##' @return xxxxxxxxxxxxxx 
##' @author Samuel Wieczorek
##' @examples
scatterplotEnrichGO_HC <- function(ego, maxRes = 10){
    
    dat <- ego@result
    nRes <- min(maxRes, nrow(dat))
    dat$GeneRatio <- unlist(lapply(dat$GeneRatio, function(x){
        as.numeric(unlist(str_split(x,'/'))[1]) / as.numeric(unlist(str_split(x,'/'))[2])
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
                        name = dat$Description)
    
    
    txt_tooltip <- NULL
    txt_tooltip <- paste("<b> p.adjust </b>: {point.pAdjust} <br> ", 
                         "<b> Count </b>: {point.z} <br> ",
                         sep="")
    
    h1 <-  highchart() %>%
        hc_chart(type = "bubble") %>%
        hc_add_series(df) %>%
        hc_legend(enabled = FALSE) %>%
        hc_xAxis(type = "category", categories = df$name)  %>%
        hc_yAxis(title = list(text = "Gene Ratio"))  %>%
        hc_tooltip(headerFormat= '',
                   pointFormat = txt_tooltip) %>%
      my_hc_ExportMenu(filename = "GOEnrich_dotplot")
    

    
    return(h1)
}