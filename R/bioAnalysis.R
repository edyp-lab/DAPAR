


getUniprotID <- function(data){
    
    id <- c("P15924","P02538", "P02768", "P08779", "Q02413", "P14923",
"P07355", "P02788", "P29508", "P63261", "Q04695", "Q8N1N4", "P01876", "P25311",
"P01833", "P06733", "Q15149", "Q01469", "P31944","Q6KB66", "P19013", "Q08188",
"Q86YZ3", "P13646", "P04083", "P02545", "P11021", "P04259", "Q06830", "P60174",
"P14618", "P31151", "P31947",  "Q96P63", "P04040", "P02787", "P04406", "P0DMV9",
"P68371", "Q8WVV4", "P13639", "P35579", "P01040", "P05089", "P01834", "P61626",
"P68363", "O75635", "P01857", "P62987", "P10599", "P00338", "O75223", "P07339",
"Q9UGM3", "Q9UI42", "O43707", "Q9NZH8", "P07900", "P01009", "P04792", "P63104",
"P00441", "P62937", "P60842", "P11142", "Q92820", "O75369", "P58107", "P47929",
"Q13867", "Q6P4A8", "P68871", "P62258", "P04075", "P00558", "O60911", "P07384",
"P18206", "Q9C075", "P31025", "P48594", "P09211", "Q6UWP8", "Q14574",  "O75342",
"Q5T750", "P0DOY2", "P20933", "P26641", "P36952", "P06396", "P25788", "P40926",
"P14735", "Q9Y6R7", "Q9HCY8", "P22528", "Q15828", "P80188", "P30740", "P01011",
"P02763", "P00738", "P18510", "A8K2U0", "P01877", "P01619", "P06870", "P31949",
"P19971", "P05090", "P40121", "P11279", "P61160", "Q9BQ50", "P23396", "P07858",
"P59998", "P23284", "P27482", "P08865", "P13473", "P47756", "P61916", "P49720",
"P42357", "P25705", "P50395", "P63244", "P48637", "Q96FX8", "P04080", "O95274",
"Q9UIV8", "P00491", "P09972", "P20930", "P19012", "P01766", "P01860", "P01764",
"P06312", "P01871") 
        return(id)
    
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





###################################
getUniprotID_FromString <- function(x){
    x <- unlist(x)
    uniprotSepIndices <- which(x=="|")
    if (length(uniprotSepIndices) == 2) { 
        x <- paste0(x, collapse="")
        res <- substr(x,uniprotSepIndices[1], uniprotSepIndices[2]-2)
        
    } else { res <- NA}
    return(res)
}

###################################
getUniprotID_FromVector <- function(dat){
    require(stringr)
    d <- str_split(dat, "|", Inf)
    uniprotID <- lapply(d,test)
    
    return(unlist(uniprotID))
}
      