

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
    
    ## ENRICHMENT : GO over-representation test
    
    gene <- bitr(data, fromType=idFrom, toType=idTo, OrgDb=orgdb)
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
# getUniprotID_FromString <- function(x){
#     x <- unlist(x)
#     uniprotSepIndices <- which(x=="|")
#     if (length(uniprotSepIndices) == 2) { 
#         x <- paste0(x, collapse="")
#         res <- substr(x,uniprotSepIndices[1], uniprotSepIndices[2]-2)
#         
#     } else { res <- NA}
#     return(res)
# }

##' xxxxxxxxxxxxxx
##' 
##' @title xxxxxxxxxxxxxx
##' @param dat xxxxxxxxxxxxxx 
##' @return xxxxxxxxxxxxxx 
##' @author Samuel Wieczorek
# getUniprotID_FromVector <- function(dat){
#     require(stringr)
#     d <- str_split(dat, "|", Inf)
#     uniprotID <- lapply(d,test)
#     
#     return(unlist(uniprotID))
# }


##' Returns the universe = totality of ID of an OrgDb annotation package
##' 
##' @title Returns the totality of ENTREZ ID (gene id) of an OrgDb annotation 
##' package
##' @param orgdb a Bioconductor OrgDb annotation package 
##' @return A vector of ENTREZ ID (totality if the ID for the package) 
##' @author Florence Combes
##' 
##' careful : org.Pf.plasmo.db : no 'ENTREZID' but 'ORF'
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
##' @param dat xxxxxxxxxxxxxx 
##' @return xxxxxxxxxxxxxx 
##' @author Samuel Wieczorek
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




barplotGroupGO_HC <- function(ggo, nRes=5){
    
    dat <- ggo@result
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


barplotEnrichGO_HC <- function(ego, nRes = 5){
    
    dat <- ego@result
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
        hc_exporting(enabled = TRUE,filename = "GOEnrich_barplot") %>%
        hc_legend(enabled = FALSE) %>%
        hc_plotOptions(bar = list(
            pointWidth=60,
            dataLabels = list(enabled = TRUE)))
    
     return(h1)
}



scatterplotEnrichGO_HC <- function(ego, nRes = 10){
    
    dat <- ego@result
    
    dat$GeneRatio <- unlist(lapply(dat$GeneRatio, function(x){
        as.numeric(unlist(str_split(x,'/'))[1]) / as.numeric(unlist(str_split(x,'/'))[2])
    }))
    
    
    n <- which(dat$Count==0)
    if (length(n) > 0){dat <- dat[-which(dat$Count==0),]}
    dat <- dat[order(dat$Count, decreasing=TRUE),]
    dat <- dat[seq(1:nRes),]
    
    
    colfunc <- colorRampPalette(c("red","royalblue"))
    nbBreaks <- 100
    series <- list()
    pal <- colfunc(nbBreaks)
    t <- log(dat$p.adjust)
    d <- (max(t) - min(t))/nbBreaks
    base <- seq(from=min(t), to=max(t), by = d)
    myColorsIndex <- unlist(lapply(t, function(x){last(which(x > base))}))
    myColorsIndex[which(is.na(myColorsIndex))] <- 1
    
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
        hc_exporting(enabled = TRUE,filename = "GOEnrich_dotplot") 
    
    
    return(h1)
}