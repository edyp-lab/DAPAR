

##' This function is a wrappper to the function groupGO from the
##' package clusterProfiler. Given a vector of genes/proteins, it returns the 
##' GO profile at a specific level. 
##' 
##' It returns a groupGOResult instance. 
##' 
##' 
##' @title Calculates the GO profile of a vector of genes/proteins at specific 
##' level
##' @param data A vector of ID (genes or proteins !!DIRE LESQUELS!!)
##' @param idFrom
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
##' GO profile at a specific level. 
##' 
##' It returns a groupGOResult instance. 
##' 
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
getUniprotID_FromString <- function(x){
    x <- unlist(x)
    uniprotSepIndices <- which(x=="|")
    if (length(uniprotSepIndices) == 2) { 
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
    uniprotID <- lapply(d,test)
    
    return(unlist(uniprotID))
}


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
    univ<-keys(get(orgdb), keytype="ENTREZID")
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


      