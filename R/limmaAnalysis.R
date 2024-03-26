#' @title Check if xxxxxx
#'
#' @param tab A data.frame which correspond to xxxxxx
#'
#' @return A list of two items
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' test.design(Biobase::pData(Exp1_R25_pept)[, seq_len(3)])
#'
#' @export
#'
test.design <- function(tab) {
  valid <- TRUE
  txt <- NULL
  level <- NULL
  
  level.a <- factor(tab[, 1], ordered = TRUE)
  level.b <- factor(tab[, 2], ordered = TRUE)
  name.level.a <- colnames(tab)[1]
  name.level.b <- colnames(tab)[2]
  
  level.c <- NULL
  if (ncol(tab) == 3) {
    level.c <- factor(tab[, 3], ordered = TRUE)
    name.level.c <- colnames(tab)[3]
  }
  
  # verification intersection sur B
  # verification de la non redondance'intersection
  # vide entre les groupes
  uniqueA <- unique(level.a)
  ll <- lapply(
    uniqueA,
    function(x) {
      as.character(level.b)[which(level.a == x)]
    }
  )
  n <- NULL
  for (i in seq_len(length(uniqueA) - 1)) {
    for (j in seq.int(from=(i + 1), to = length(uniqueA))) {
      n <- c(n, intersect(ll[[i]], ll[[j]]))
    }
  }
  
  if (length(n) > 0) {
    valid <- FALSE
    txt <- c(txt, paste0(
      "The value ",
      n,
      " in column '",
      colnames(tab)[2],
      "' is not correctly set.\n"
    ))
  }
  
  
  # verification si niveau hierarchique inf
  if (length(levels(level.a)) == length(levels(level.b))) {
    ## c'est un design de niveau n-1 en fait
    valid <- FALSE
    txt <- c(
      txt,
      paste0(
        "The column ",
        name.level.b,
        " is not informative. ",
        "Thus, the design is not of level (n-1).\n"
      )
    )
  } else if (!is.null(level.c)) {
    if (length(levels(level.b)) == length(levels(level.c))) {
      ## c'est un design de niveau n-1 en fait
      valid <- FALSE
      txt <- c(
        txt,
        paste0(
          "The column ",
          name.level.c,
          " is not informative. ",
          "Thus, the design is of level (n-1).\n"
        )
      )
    }
  }
  
  # verification si niveau non informatif
  return(list(
    valid = valid,
    warn = txt
  ))
}




#' @title Check if the design is valid
#'
#' @param conds A vector
#'
#' @return A list
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' check.conditions(Biobase::pData(Exp1_R25_pept)$Condition)
#'
#' @export
#'
check.conditions <- function(conds) {
  res <- list(valid = TRUE, warn = NULL)
  
  if (("" %in% conds) || (NA %in% conds)) {
    res <- list(valid = FALSE, 
                warn = "The conditions are note full filled.")
    return(res)
  }
  
  # Check if there is at least two conditions
  if (length(unique(conds)) < 2) {
    res <- list(valid = FALSE, 
                warn = "The design must contain at least two conditions.")
    return(res)
  }
  
  
  # check if each condition has at least two values
  nValPerCond <- unlist(lapply(unique(conds), function(x) {
    length(conds[which(conds == x)])
  }))
  if (all(nValPerCond < 2)) {
    res <- list(valid = FALSE, 
                warn = "The design must contain at least two values per condition.")
    return(res)
  }
  
  return(res)
}



#' @title Check if the design is valid
#'
#' @param sTab The data.frame which correspond to the `pData()` function 
#' of package `MSnbase`.
#'
#' @return A boolean
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' check.design(Biobase::pData(Exp1_R25_pept)[, seq_len(3)])
#'
#' @export
#'
check.design <- function(sTab) {
  res <- list(valid = FALSE, warn = NULL)
  
  names <- colnames(sTab)
  level.design <- ncol(sTab) - 2
  
  
  res <- check.conditions(sTab$Condition)
  if (!res$valid) {
    return(res)
  }
  # Check if all the column are fullfilled
  
  if (level.design == 1) {
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Bio.Rep colmumn are not full filled.")
      return(res)
    }
  } else if (level.design == 2) {
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Bio.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Tech.Rep) || (NA %in% sTab$Tech.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Tech.Rep colmumn are not full filled.")
      return(res)
    }
  } else if (level.design == 3) {
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Bio.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Tech.Rep) || (NA %in% sTab$Tech.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Tech.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Analyt.Rep) || (NA %in% sTab$Analyt.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Analyt.Rep colmumn are not full filled.")
      return(res)
    }
  }
  
  # Check if the hierarchy of the design is correct
  if (level.design == 1) {
    res <- test.design(sTab[, c("Condition", "Bio.Rep")])
  } else if (level.design == 2) {
    res <- test.design(sTab[, c("Condition", "Bio.Rep", "Tech.Rep")])
  } else if (level.design == 3) {
    res <- test.design(sTab[, c("Condition", "Bio.Rep", "Tech.Rep")])
    if (res$valid) {
      res <- test.design(sTab[, c("Bio.Rep", "Tech.Rep", "Analyt.Rep")])
    }
  }
  
  return(res)
}





#' @title Builds the design matrix
#'
#' @param sTab The data.frame which correspond to the `pData()` function 
#' of package `MSnbase`.
#'
#' @return A design matrix
#'
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' make.design(Biobase::pData(Exp1_R25_pept))
#'
#' @export
#' 
make.design <- function(sTab) {
  if (!check.design(sTab)$valid) {
    warning("The design matrix is not correct.")
    warning(check.design(sTab)$warn)
    return(NULL)
  }
  
  n <- ncol(sTab)
  if (n == 1 || n == 2) {
    stop.txt <- paste0("Error in design matrix dimensions which must ", 
                       "have at least 3 columns.")
    stop(stop.txt)
  }
  
  res <- do.call(paste0("make.design.", (n - 2)), list(sTab))
  
  return(res)
}



#' @title Builds the design matrix for designs of level 1
#'
#' @param sTab The data.frame which correspond to the `pData()` function 
#' of package `MSnbase`.
#'
#' @return A design matrix
#'
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' make.design.1(Biobase::pData(Exp1_R25_pept))
#'
#' @export
make.design.1 <- function(sTab) {
  Conditions <- factor(sTab$Condition, ordered = TRUE)
  nb_cond <- length(unique(Conditions))
  nb_samples <- nrow(sTab)
  
  # CGet the number of replicates per condition
  nb_Rep <- rep(0, nb_cond)
  for (i in seq_len(nb_cond)) {
    nb_Rep[i] <- sum((Conditions == unique(Conditions)[i]))
  }
  
  design <- matrix(0, nb_samples, nb_cond)
  n0 <- 1
  coln <- NULL
  
  for (j in seq_len(nb_cond)) {
    # When design == 1, one uses letters to have up to 26 possible conditions
    if (getDesignLevel(sTab)==1)
      id <- LETTERS[j]
    else
      id <- j
    
    coln <- c(coln, paste("Condition", id, collapse = NULL, sep = ""))
    seq <- seq.int(from=n0, to=(n0 + nb_Rep[j] - 1))
    design[seq, j] <- rep(1, length(seq))
    n0 <- n0 + nb_Rep[j]
  }
  
  
  
  
  #for (j in seq_len(nb_cond)) {
  #  coln <- c(coln, paste("Condition", j, collapse = NULL, sep = ""))
  #  design[seq.int(from=n0, to=(n0 + nb_Rep[j] - 1)), j] <- rep(1, 
  #                                                             length(seq.int(from=n0, to = (n0 + nb_Rep[j] - 1)))
  # )
  #  n0 <- n0 + nb_Rep[j]
  # }
  colnames(design) <- coln
  
  return(design)
}





#' @title Builds the design matrix for designs of level 2
#'
#' @param sTab The data.frame which correspond to the `pData()` function 
#' of package `MSnbase`.
#'
#' @return A design matrix
#'
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package='DAPARdata')
#' make.design.2(Biobase::pData(Exp1_R25_pept))
#'
#'
#' @export
#'
make.design.2 <- function(sTab) {
  pkgs.require('stats')
  
  
  Condition <- factor(sTab$Condition, levels = unique(sTab$Condition))
  RepBio <- factor(sTab$Bio.Rep, levels = unique(sTab$Bio.Rep))
  
  # Renome the levels of factor
  levels(Condition) <- seq_len(length(levels(Condition)))
  levels(RepBio) <- seq_len(length(levels(RepBio)))
  
  # Initial design matrix
  df <- rep(0, nrow(sTab))
  names(df) <- rownames(sTab)
  design <- stats::model.matrix(df ~ 0 + Condition:RepBio)
  
  # Remove empty columns in the design matrix
  design <- design[, (apply(design, 2, sum) > 0)]
  # Remove identical columns in the design matrix
  coldel <- -1
  for (i in seq_len(length(design[1, ]) - 1)) {
    d2 <- as.matrix(
      design[, seq.int(from=(i + 1), to = length(design[1, ]))]
    )
    for (j in seq_len(length(d2[1, ]))) {
      d2[, j] <- d2[, j] - design[, i]
    }
    e <- as.matrix(stats::rnorm(length(design[, 1]), 10, 1))
    sd2 <- t(e) %*% d2
    liste <- which(sd2 == 0)
    coldel <- c(coldel, liste + i)
  }
  design <- design[, (seq_len(length(design[1, ]))) != coldel]
  colnames(design) <- make.names(colnames(design))
  return(design)
}





#' @title Builds the design matrix for designs of level 3
#'
#' @param sTab The data.frame which correspond to the `pData()` function 
#' of package `MSnbase`.
#'
#' @return A design matrix
#'
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' sTab <- cbind(Biobase::pData(Exp1_R25_pept), Tech.Rep = 1:6)
#' make.design.3(sTab)
#'
#'
#' @export
#'
make.design.3 <- function(sTab) {
  
  pkgs.require('stats')
  
  Condition <- factor(sTab$Condition, levels = unique(sTab$Condition))
  RepBio <- factor(sTab$Bio.Rep, levels = unique(sTab$Bio.Rep))
  RepTech <- factor(sTab$Tech.Rep, levels = unique(sTab$Tech.Rep))
  
  
  # Rename the levels of factor
  levels(Condition) <- seq_len(length(levels(Condition)))
  levels(RepBio) <- seq_len(length(levels(RepBio)))
  levels(RepTech) <- seq_len(length(levels(RepTech)))
  
  
  # Initial design matrix
  df <- rep(0, nrow(sTab))
  names(df) <- rownames(sTab)
  design <- stats::model.matrix(df ~ 0 + Condition:RepBio:RepTech)
  
  # Remove empty columns in the design matrix
  design <- design[, (apply(design, 2, sum) > 0)]
  
  # Remove identical columns in the design matrix
  coldel <- -1
  for (i in seq_len(length(design[1, ]) - 1)) {
    d2 <- as.matrix(design[, seq.int(from = (i + 1), 
                                     to = length(design[1, ]))])
    for (j in seq_len(length(d2[1, ]))) {
      d2[, j] <- d2[, j] - design[, i]
    }
    e <- as.matrix(stats::rnorm(length(design[, 1]), 10, 1))
    sd2 <- t(e) %*% d2
    liste <- which(sd2 == 0)
    coldel <- c(coldel, liste + i)
  }
  design <- design[, seq_len(length(design[1, ])) != coldel]
  colnames(design) <- make.names(colnames(design))
  return(design)
}


#' @title xxx
#' @description xxx
#' 
#' @param sTab xxx
#' 
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' sTab <- Biobase::pData(Exp1_R25_pept)
#' getDesignLevel(sTab)
#'
#' @export
#' 
getDesignLevel <- function(sTab){
  
  level <- ncol(sTab) - 2
  
  return (level)
}




#' @title Builds the contrast matrix
#'
#' @param design The data.frame which correspond to the `pData()` function 
#' of package `MSnbase`.
#'
#' @param condition xxxxx
#'
#' @param contrast An integer that Indicates if the test consists of the
#' comparison of each biological condition versus each of the other ones
#' (Contrast=1; for example H0:"C1=C2" vs H1:"C1!=C2", etc.)
#' or each condition versus all others (Contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
#'  H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
#' @param design.level xxx
#'
#' @return A constrat matrix
#'
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package='DAPARdata')
#' design <- make.design(Biobase::pData(Exp1_R25_pept))
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' make.contrast(design, conds)
#'
#' @export
#'
make.contrast <- function(design, 
                          condition, 
                          contrast = 1,
                          design.level = 1) {
  
  
  aggreg.column.design <- function(design, 
                                   Condition,
                                   design.level) {
    nb.cond <- length(unique(Condition))
    name.col <- colnames(design)
    name.cond <- NULL
    nb.col <- NULL
    for (i in seq_len(nb.cond)) {
      col.select <- NULL
      id <- if (design.level == 1) LETTERS[i] else i
      
      col.name.begin <- paste("Condition", id, sep = "")
      nc <- nchar(col.name.begin)
      for (j in seq_len(length(design[1, ]))) {
        if (substr(name.col[j], 1, nc) == col.name.begin) {
          col.select <- c(col.select, j)
        }
      }
      name.aggreg <- NULL
      for (j in seq_len(length(col.select))) {
        name.aggreg <- paste(name.aggreg, 
                             name.col[col.select[j]], 
                             sep = "+")
      }
      name.aggreg <- substr(name.aggreg, 2, nchar(name.aggreg))
      name.cond <- c(name.cond, name.aggreg)
      nb.col <- c(nb.col, length(col.select))
    }
    return(list(name.cond, nb.col))
  }
  
  
  
  nb.cond <- length(unique(condition))
  r <- aggreg.column.design(design, condition, design.level)
  label.agg <- r[[1]]
  nb.agg <- r[[2]]
  k <- 1
  
  if (contrast == 1) {
    ## Contrast for One vs One
    contra <- rep(0, sum(seq_len((nb.cond - 1))))
    for (i in seq_len(nb.cond - 1)) {
      for (j in seq.int(from = (i + 1), to = nb.cond)) {
        contra[k] <- c(paste(
          "(", label.agg[i], ")/",
          nb.agg[i], "-(", label.agg[j], ")/",
          nb.agg[j]
        ))
        k <- k + 1
      }
    }
  } else if (contrast == 2) {
    ## Contrast for One vs All
    contra <- rep(0, nb.cond)
    for (i in seq_len(nb.cond)) {
      contra[k] <- c(paste("(", label.agg[i], ")/", nb.agg[i]))
      nb <- sum(nb.agg[seq_len(nb.cond)[seq_len(nb.cond) != i]])
      for (j in seq_len(nb.cond)[seq_len(nb.cond) != i]) {
        contra[k] <- c(paste(contra[k], "-(", label.agg[j], ")/", nb))
      }
      k <- k + 1
    }
  }
  
  return(contra)
}


#' @title Computes a hierarchical differential analysis
#'
#' @param qData A matrix of quantitative data, without any missing values.
#'
#' @param sTab A dataframe of experimental design (Biobase::pData()).
#'
#' @param comp.type A string that corresponds to the type of comparison.
#' Values are: 'anova1way', 'OnevsOne' and 'OnevsAll'; default is 'OnevsOne'.
#'
#' @return A list of two dataframes : logFC and P_Value. The first one contains
#' the logFC values of all the comparisons (one column for one comparison),
#' the second one contains the pvalue of all the comparisons (one column for
#' one comparison). The names of the columns for those two dataframes
#' are identical and correspond to the description of the comparison.
#'
#' @author Hélène Borges, Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept
#' qData <- Biobase::exprs(obj)
#' sTab <- Biobase::pData(obj)
#' limma <- limmaCompleteTest(qData, sTab, comp.type = "anova1way")
#'
#' @export
#'
#'
limmaCompleteTest <- function(qData, sTab, comp.type = "OnevsOne") {
  
  pkgs.require(c('dplyr', 'limma', 'tidyr'))
  
  
  level <- getDesignLevel(sTab)
  if (level == 1 && length(unique(sTab$Condition)) > 26){
    message('DAPAR with Limma does not handle datasets with more 26 conditions 
            in a 1-level design')
    return(NULL)
  } else if ((level == 2 || level == 3) && length(unique(sTab$Condition)) >= 10){
    message('DAPAR with Limma does not handle datasets with more 9 conditions 
            in a 2-level and 3-level designs')
    return(NULL)
  }
  
  switch(comp.type,
         OnevsOne = contrast <- 1,
         OnevsAll = contrast <- 2
  )
  #sTab.old <- sTab
  conds <- factor(sTab$Condition, levels = unique(sTab$Condition))
  #sTab <- sTab[unlist(lapply(split(sTab, conds), function(x) {x["Sample.name"]})), ]
  #qData <- qData[, unlist(lapply(split(sTab.old, conds), function(x) {x["Sample.name"]}))]
  #conds <- conds[order(conds)]
  
  res.l <- NULL
  
  design.matrix <- make.design(sTab)
  
  if (!is.null(design.matrix)) {
    if (comp.type == "OnevsOne" || comp.type == "OnevsAll") {
      contra <- make.contrast(design.matrix, 
                              condition = conds, 
                              contrast,
                              design.level = getDesignLevel(sTab))
      cmtx <- limma::makeContrasts(
        contrasts = contra,
        levels = make.names(colnames(design.matrix))
      )
      fit <- limma::eBayes(
        limma::contrasts.fit(limma::lmFit(qData, design.matrix), cmtx))
      res.l <- formatLimmaResult(fit, 
                                 conds, 
                                 contrast,
                                 design.level = getDesignLevel(sTab))
    } else if (comp.type == "anova1way") {
      # make the orthogonal contrasts
      contrasts <- tidyr::crossing(
        A = colnames(design.matrix), 
        B = colnames(design.matrix), 
        .name_repair = "minimal") %>%
        dplyr::filter(A != B)
      orthogonal_contrasts <- dplyr::filter(
        contrasts, !duplicated(paste0(pmax(A, B), pmin(A, B)))) %>%
        dplyr::mutate(contrasts = stringr::str_glue("{A}-{B}"))
      # make the contrasts in a format adapted for limma functions
      contrasts_limma_format <- limma::makeContrasts(
        contrasts = orthogonal_contrasts$contrasts,
        levels = colnames(design.matrix)
      )
      ebayes_fit <- limma::eBayes(
        limma::contrasts.fit(
          limma::lmFit(qData, design.matrix), 
          contrasts_limma_format)
      )
      fit_table <- limma::topTable(
        ebayes_fit, 
        sort.by = "none", 
        number = nrow(qData)
      )
      fit_pvalue <- dplyr::select(fit_table, "anova_1way_pval" = P.Value)
      res.l <- list(
        "logFC" = data.frame(
          "anova_1way_logFC" = matrix(NA, nrow = nrow(fit_pvalue)),
          row.names = rownames(fit_pvalue)
        ),
        "P_Value" = fit_pvalue
      )
    }
  }
  return(res.l)
}





#' @title xxxx
#'
#' @param fit xxxx
#'
#' @param conds xxxx
#'
#' @param contrast xxxx
#' @param design.level xxx
#'
#' @return A list of two dataframes : logFC and P_Value. The first one contains
#' the logFC values of all the comparisons (one column for one comparison),
#' the second one contains the pvalue of all the comparisons (one column for
#' one comparison). The names of the columns for those two dataframes are
#' identical and correspond to the description of the comparison.
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
#' limma <- limmaCompleteTest(qData, sTab)
#'
formatLimmaResult <- function(fit, conds, contrast, design.level) {
  pkgs.require('stringr')
  
  
  res <- cbind(fit$coefficients, fit$p.value)
  
  # how many comparisons have been done (and thus how many columns of pval)
  Compa.Nb <- dim(fit$p.value)[2]
  # empty colnames vector
  cn <- c()
  for (i in seq_len(Compa.Nb)) {
    
    # not the same syntax to pars if Contast=1 or Contrast=2
    if (contrast == 1) {
      if(design.level == 1) {
        compa <- stringr::str_match_all(
          colnames(fit$p.value)[i], 
          "[[:space:]]Condition([[:letter:]]+)")[[1]]
        
        tmp1 <- unique(conds)[which(LETTERS == compa[1, 2])]
        tmp2 <- unique(conds)[which(LETTERS == compa[2, 2])]
        cn[i] <- paste(tmp1, "_vs_", tmp2, sep = "")
        
      } else {
        compa <- stringr::str_match_all(
          colnames(fit$p.value)[i], 
          "[[:space:]]Condition([[:digit:]]+)")[[1]]
        
        tmp1 <- unique(conds)[as.numeric(compa[1, 2])]
        tmp2 <- unique(conds)[as.numeric(compa[2, 2])]
        cn[i] <- paste(tmp1, "_vs_", tmp2, sep = "")
      }
    }
    if (contrast == 2) {
      
      if(design.level == 1) {
        compa <- stringr::str_match_all(
          colnames(fit$p.value)[i], 
          "[[:space:]]Condition([[:letter:]]+)")[[1]]
        
        tmp <- unique(conds)[which(LETTERS == compa[1, 2])]
        cn[i] <- paste(tmp, "_vs_(all-", tmp, ")", sep = "")
        
      } else {
        compa <- stringr::str_match_all(
          colnames(fit$p.value)[i], 
          "[[:space:]]Condition([[:digit:]]+)")[[1]]
        
        tmp <- unique(conds)[as.numeric(compa[1, 2])]
        cn[i] <- paste(tmp, "_vs_(all-", tmp, ")", sep = "")
      }
      
      
      
    }
  }
  
  
  res.l <- list(
    logFC = as.data.frame(res[, seq_len(Compa.Nb)]),
    P_Value = as.data.frame(res[, -(seq_len(Compa.Nb))])
  )
  
  colnames(res.l$logFC) <- gsub("[ ]", "", paste(cn, "logFC", sep = "_"))
  colnames(res.l$P_Value) <- gsub("[ ]", "", paste(cn, "pval", sep = "_"))
  
  return(res.l)
}
