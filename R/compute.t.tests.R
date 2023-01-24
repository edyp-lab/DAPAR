
#' @title xxxxxx
#'
#' @param obj A matrix of quantitative data, without any missing values.
#'
#' @param contrast Indicates if the test consists of the comparison of each
#' biological condition versus
#' each of the other ones (contrast=1;
#' for example H0:"C1=C2" vs H1:"C1!=C2", etc.)
#' or each condition versus all others (contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
#' H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
#'
#' @param type xxxxx
#'
#' @return A list of two items : logFC and P_Value; both are dataframe. The
#' first one contains the logFC values of all the comparisons (one column for
#' one comparison), the second one contains the pvalue of all the comparisons
#' (one column for one comparison). The names of the columns for those two
#' dataframes are identical and correspond to the description of the comparison.
#'
#' @author Florence Combes, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(1000)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(GetMetacell(obj), "Missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' ttest <- compute_t_tests(obj$new)
#'
#'
#' @export
#'
compute_t_tests <- function(obj, contrast = "OnevsOne", type = "Student") {
  
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("Please install stats: BiocManager::install('stats')")
  }
  
  
  switch(type,
         Student = .type <- TRUE,
         Welch = .type <- FALSE
  )
  
  
  qData <- Biobase::exprs(obj)
  sTab <- Biobase::pData(obj)
  res <- list()
  logFC <- list()
  P_Value <- list()
  
  nbComp <- NULL
  
  #sTab.old <- sTab
  Conditions.f <- factor(sTab$Condition, levels = unique(sTab$Condition))
  #sTab <- sTab[unlist(lapply(split(sTab, Conditions.f), function(x) {x["Sample.name"]})), ]
  #qData <- qData[, unlist(lapply(split(sTab.old, Conditions.f), function(x) { x["Sample.name"]}))]
  Conditions <- sTab$Condition
  
  
  # Cond<-levels(Conditions.f)
  Cond.Nb <- length(levels(Conditions.f))
  
  
  if (contrast == "OnevsOne") {
    res.tmp <- NULL
    nbComp <- Cond.Nb * (Cond.Nb - 1) / 2
    
    for (i in seq_len(Cond.Nb - 1)) {
      for (j in seq.int(from = (i + 1), to = Cond.Nb)) {
        
        c1Indice <- which(Conditions == levels(Conditions.f)[i])
        c2Indice <- which(Conditions == levels(Conditions.f)[j])
        
        res.tmp <- apply(
          qData[, c(c1Indice, c2Indice)], 1,
          function(x) {
            stats::t.test(x ~ Conditions[c(c1Indice, c2Indice)], var.equal = .type)
          }
        )
        p.tmp <- unlist(lapply(res.tmp, function(x) x$p.value))
        m1.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(x$estimate[1])))
        m2.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(x$estimate[2])))
        m1.name <- names(unlist(lapply(res.tmp, function(x) x$estimate[1])))[1]
        m2.name <- names(unlist(lapply(res.tmp, function(x) x$estimate[2])))[1]
        logFC.tmp <- m1.tmp - m2.tmp
        
        #if (grepl(levels(Conditions.f)[i], m2.name)) {
        #     logFC.tmp <- -logFC.tmp
        # }
        
        txt <- paste(levels(Conditions.f)[i], "_vs_", 
                     levels(Conditions.f)[j], sep = "")
        
        logFC[[paste(txt, "logFC", sep = "_")]] <- logFC.tmp
        P_Value[[paste(txt, "pval", sep = "_")]] <- p.tmp
      }
    }
  } ## end Contrast==1
  
  if (contrast == "OnevsAll") {
    nbComp <- Cond.Nb
    
    for (i in seq_len(nbComp)) {
      c1 <- which(Conditions == levels(Conditions.f)[i])
      
      Cond.t.all <- seq_len(length(Conditions))
      Cond.t.all[c1] <- levels(Conditions.f)[i]
      Cond.t.all[-c1] <- "all"
      
      res.tmp <- apply(qData, 1,function(x) {
        stats::t.test(x ~ Cond.t.all, var.equal = .type)
      }
      )
      
      p.tmp <- unlist(lapply(res.tmp, function(x) x$p.value))
      m1.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(x$estimate[1])))
      m2.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(x$estimate[2])))
      m1.name <- names(unlist(lapply(res.tmp, function(x) x$estimate[1])))[1]
      m2.name <- names(unlist(lapply(res.tmp, function(x) x$estimate[2])))[1]
      logFC.tmp <- m1.tmp - m2.tmp
      #if (grepl(levels(Conditions.f)[i], m2.name)) {
      #    logFC.tmp <- -logFC.tmp
      #}
      
      txt <- paste(levels(Conditions.f)[i], "_vs_(all-", 
                   levels(Conditions.f)[i], ")", sep = "")
      
      logFC[[paste(txt, "logFC", sep = "_")]] <- logFC.tmp
      P_Value[[paste(txt, "pval", sep = "_")]] <- p.tmp
    }
  } # End Contrast=2
  
  
  res.l <- list(
    logFC = as.data.frame(logFC),
    P_Value = as.data.frame(P_Value)
  )
  colnames(res.l$logFC) <- names(logFC)
  colnames(res.l$P_Value) <- names(P_Value)
  
  return(res.l)
}
