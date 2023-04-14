#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param current_protein xxx
#' @param conditions xxx
#'
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
OWAnova <- function(current_protein, conditions){
  intensities <- unname(unlist(current_protein))
  aov_model <- aov(formula = intensities ~ conditions, data = NULL)
  return(aov_model)
}

#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param obj xxx
#' '
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
applyAnovasOnProteins <- function(obj){
  qData <- Biobase::exprs(obj)
  sTab <- Biobase::pData(obj)
  anova_models <- t(apply(qData,1, OWAnova, conditions=as.factor(sTab$Condition)))
  names(anova_models) <- rownames(qData)
  return(anova_models)
}

#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param aov_fits xxx
#' @param test xxx
#'
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
testAnovaModels <- function(aov_fits, test = "Omnibus"){  
  if (! requireNamespace("multcomp", quietly = TRUE)) {
    stop("Please install multcomp: BiocManager::install('multcomp')")
  }
  switch(test,
         Omnibus={
           omnibus_tests_summaries <- t(sapply(aov_fits,
                                               function(x) unlist(summary(x))
           ))
           res <- list("logFC" = data.frame("anova_1way_logFC" = matrix(NA, nrow = length(aov_fits)), row.names = names(aov_fits)),
                       "P_Value" = data.frame("anova_1way_pval" = omnibus_tests_summaries[,9], row.names = names(aov_fits)))
         }, TukeyHSD={
           tukeyHSD_tests_summaries <- lapply(aov_fits, 
                                              function(x) TukeyHSD(x, which = "conditions")$conditions)
           res <- formatHSDResults(tukeyHSD_tests_summaries)
         }, TukeySinglestep={
           tukeySS_tests_summaries <- lapply(aov_fits, 
                                             function(x) summary(multcomp::glht(x, linfct = multcomp::mcp(conditions = "Tukey")),
                                                                 test = multcomp::adjusted("single-step"))
           )
           res <- formatPHTResults(tukeySS_tests_summaries)
         }, TukeyStepwise={
           tukeySW_tests_summaries <- lapply(aov_fits, 
                                             function(x) summary(multcomp::glht(x, linfct = multcomp::mcp(conditions = "Tukey")),
                                                                 test = multcomp::adjusted("Westfall"))
           )
           res <- formatPHTResults(tukeySW_tests_summaries)
         }, TukeyNoMTC={
           dunnettSW_tests_summaries <- lapply(aov_fits,
                                               function(x) summary(multcomp::glht(x, linfct = multcomp::mcp(conditions = "Tukey")),
                                                                   test = multcomp::adjusted("none"))
           )
           res <- formatPHTResults(dunnettSW_tests_summaries)
         }, DunnettSinglestep={
           dunnettSS_tests_summaries <- lapply(aov_fits,
                                               function(x) summary(multcomp::glht(x, linfct = multcomp::mcp(conditions = "Dunnett")),
                                                                   test = multcomp::adjusted("single-step"))
           )
           res <- formatPHTResults(dunnettSS_tests_summaries)
         }, DunnettStepwise={
           dunnettSW_tests_summaries <- lapply(aov_fits,
                                               function(x) summary(multcomp::glht(x, linfct = multcomp::mcp(conditions = "Dunnett")),
                                                                   test = multcomp::adjusted("free"))
           )
           res <- formatPHTResults(dunnettSW_tests_summaries)
         }, DunnettNoMTC={
           dunnettSW_tests_summaries <- lapply(aov_fits,
                                               function(x) summary(multcomp::glht(x, linfct = multcomp::mcp(conditions = "Dunnett")),
                                                                   test = multcomp::adjusted("none"))
           )
           res <- formatPHTResults(dunnettSW_tests_summaries)
         }
  )
  return(res)
}

#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param post_hoc_models_summaries xxx
#'
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
formatHSDResults <- function(post_hoc_models_summaries){
  if (! requireNamespace("purrr", quietly = TRUE)) {
    stop("Please install purrr: BiocManager::install('purrr')")
  }
  # get the fold-changes
  res_coeffs <- lapply(post_hoc_models_summaries, function(x) x[,1])
  logFC <- data.frame(purrr::map_dfr(res_coeffs, cbind),
                      row.names = names(post_hoc_models_summaries[[1]][,1]))
  logFC <- as.data.frame(t(logFC))
  # extract the FWER ajusted p-values
  res_pvals <- lapply(post_hoc_models_summaries, function(x) x[,4])
  pvals <- data.frame(purrr::map_dfr(res_pvals, cbind),
                      row.names = names(post_hoc_models_summaries[[1]][,4]))
  pvals <- as.data.frame(t(pvals))
  res <- list("logFC" = logFC,
              "P_Value" = pvals)
  # formatting of column names for consistency with the limma and t-test code
  colnames(res$logFC) <- stringr::str_replace(colnames(res$logFC), "-", "_vs_")
  colnames(res$P_Value) <- stringr::str_replace(colnames(res$P_Value), "-", "_vs_")
  colnames(res$logFC) <- stringr::str_c(colnames(res$logFC), "_logFC")
  colnames(res$P_Value) <- stringr::str_c(colnames(res$P_Value), "_pval")
  return(res)
}

#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param post_hoc_models_summaries xxx
#'
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
formatPHTResults <- function(post_hoc_models_summaries){
  if (! requireNamespace("purrr", quietly = TRUE)) {
    stop("Please install purrr: BiocManager::install('purrr')")
  }
  # récupérer les différences entre les moyennes
  res_coeffs <- lapply(post_hoc_models_summaries, function(x) x$test$coefficients)
  logFC <- data.frame(purrr::map_dfr(res_coeffs, cbind),
                      row.names = names(post_hoc_models_summaries[[1]]$test$coefficients))
  logFC <- as.data.frame(t(logFC))
  # extract raw p-values (non-adjusted)
  res_pvals <- lapply(post_hoc_models_summaries, function(x) x$test$pvalues)
  pvals <- data.frame(purrr::map_dfr(res_pvals, cbind),
                      row.names = names(post_hoc_models_summaries[[1]]$test$coefficients))
  pvals <- as.data.frame(t(pvals))
  res <- list("logFC" = logFC,
              "P_Value" = pvals)
  # formatting of column names for consistency with the limma and t-test code
  colnames(res$logFC) <- stringr::str_replace(colnames(res$logFC), " - ", "_vs_")
  colnames(res$P_Value) <- stringr::str_replace(colnames(res$P_Value), " - ", "_vs_")
  colnames(res$logFC) <- stringr::str_c(colnames(res$logFC), "_logFC")
  colnames(res$P_Value) <- stringr::str_c(colnames(res$P_Value), "_pval")
  return(res)
}

#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param x xxx
#' @param pval.T xxx
#' @param M xxx
#'
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
thresholdpval4fdr <- function(x, pval.T, M){
  index <- which(x< pval.T)
  R <- length(index)/length(x)
  print(R)
  res <- rep(1, length(x))
  res[index] <- (adjust.p(x[index], pi0.method = M)$adjp$adjusted.p)
  return(res)
}

#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param x xxx
#' @param pval.threshold xxx
#' @param method xxx
#'
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
separateAdjPval <- function(x, pval.threshold=1.05, method=1){
  if(pval.threshold > 1){
    res <- apply(x, 2, function(x)adjust.p(x, pi0.method = method)$adjp$adjusted.p)
  } else{
    res <- as.data.frame(apply(x, 2, function(x)thresholdpval4fdr(x, pval.T = pval.threshold,  M=method)))
  }
  return(res)
}

#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param x xxx
#' @param pval.threshold xxx
#' @param method xxx
#' @param display xxx
#'
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
globalAdjPval <- function(x, pval.threshold=1.05, method=1, display = T){
  res <- x
  vec <- stack(x)$values
  index <- which(vec< pval.threshold)
  if(display)calibration.plot(vec[index], pi0.method="ALL")
  vec[index] <- adjust.p(vec[index], pi0.method = method)$adjp$adjusted.p
  vec[-index] <- 1
  res[,] <- vec
  return(res)
}

#' @title xxx
#'
#' @author Thomas Burger
#'
#' @param x xxx
#'
#' @param fdr.threshold xxx
#'
#' @return xxx
#'
#' @examples
#' NULL
#'
#' @export
#' 
compute.selection.table <- function(x, fdr.threshold){
  selection.table <- 1*(x<fdr.threshold)
  tt <- sum(selection.table)
  names(tt) <- "Total Ctr"
  pc <- sum(as.logical(rowSums(selection.table)))
  #names(pc) <- paste("DA_prot/", dim(selection.table)[1], sep="")
  names(pc) <- "DA prot"
  res <- c(colSums(selection.table), tt, pc)
}







#' @title Function to perform a One-way Anova statistical test on a MsnBase 
#' dataset
#'
#' @author Hélène Borges
#'
#' @param current_line The line currently treated from the quantitative data to
#'  perform the ANOVA
#'
#' @param conditions The conditions represent the different classes of the
#' studied factor
#'
#' @return A named vector containing all the different values of the aov model
#'
#' @examples
#' \dontrun{examples/ex_classic1wayAnova.R}
#'
#' @export
#'
classic1wayAnova <- function(current_line, conditions) {
  .Deprecated("OWAnova")
  pkgs.require('stats')
  
  
  # vector containing the protein/peptide intensities
  intensities <- unname(unlist(current_line))
  # anova on these two vectors
  aov_test <- stats::aov(formula = intensities ~ conditions, data = NULL)
  return(aov_test)
}

#' @title Wrapper for One-way Anova statistical test
#'
#' @author Hélène Borges
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param with_post_hoc a character string with 2 possible values: "Yes" and
#' "No" (default) saying if function must perform a Post-Hoc test or not.
#'
#' @param post_hoc_test character string, possible values are "No" (for no
#' test; default value) or TukeyHSD" or "Dunnett". See details of
#' \code{postHocTest()} function to choose the appropriate one.
#'
#' @details This function allows to perform a 1-way Analysis of Variance. Also
#' computes the post-hoc tests if the \code{with_post_hoc} parameter is set to
#' yes. There are two possible post-hoc tests: the Tukey Honest Significant
#' Differences (specified as "TukeyHSD") and the Dunnett test
#' (specified as "Dunnett").
#'
#' @return A list of two dataframes. First one called "logFC" contains
#' all pairwise comparisons logFC values (one column for one comparison) for
#' each analysed feature (Except in the case without post-hoc testing, for
#' which NAs are returned.); The second one named "P_Value" contains the
#' corresponding p-values.
#'
#' @examples 
#' \dontrun{examples/ex_wrapperClassic1wayAnova.R}
#'
#' @seealso [postHocTest()]
#'
#' @export
#'
wrapperClassic1wayAnova <- function(obj, 
                                    with_post_hoc = "No", 
                                    post_hoc_test = "No") {
  
  .Deprecated("testAnovaModels")
  pkgs.require('dplyr')
  
  
  qData <- Biobase::exprs(obj)
  sTab <- Biobase::pData(obj)
  if (with_post_hoc == "No") {
    anova_tests <- as.data.frame(
      t(
        apply(
          qData, 1, 
          function(x) 
            unlist(
              summary(
                classic1wayAnova(x, 
                                 conditions = as.factor(sTab$Condition)
                )
              )
            )
        )
      )
    )
    results <- dplyr::select(anova_tests, `Pr(>F)1`)
    to_return <- list(
      "logFC" = data.frame(
        "anova_1way_logFC" = matrix(NA, nrow = nrow(results)), 
        row.names = rownames(results)
      ),
      "P_Value" = data.frame(
        "anova_1way_pval" = results['Pr(>F)1'], 
        row.names = rownames(results)
      )
    )
  } else if (with_post_hoc == "Yes") {
    if (post_hoc_test == "No") {
      stop("You want to perform a post-hoc test but did not specify 
                which test. Please choose between Dunnett or TukeyHSD")
    } else if (post_hoc_test == "TukeyHSD") {
      anova_tests <- t(
        apply(qData, 1, classic1wayAnova, 
              conditions = as.factor(sTab$Condition))
      )
      names(anova_tests) <- rownames(qData)
      to_return <- postHocTest(
        aov_fits = anova_tests, 
        post_hoc_test = post_hoc_test
      )
    } else if (post_hoc_test == "Dunnett") {
      anova_tests <- t(
        apply(qData, 1, classic1wayAnova, 
              conditions = as.factor(sTab$Condition))
      )
      names(anova_tests) <- rownames(qData)
      to_return <- postHocTest(
        aov_fits = anova_tests, 
        post_hoc_test = post_hoc_test
      )
    }
  } else {
    stop("Wrong with_post_hoc parameter. Please choose between No or Yes.")
  }
  
  return(to_return)
}



#' Extract logFC and raw pvalues from multiple post-hoc models summaries
#'
#' @param post_hoc_models_summaries a list of summaries of post-hoc models.
#'
#' @return a list of 2 dataframes containing the logFC values and pvalues for
#' each comparison.
#'
#' @author Hélène Borges
#'
#' @examples
#' \dontrun{examples/ex_formatPHResults.R}
#'
#' @export
#'
formatPHResults <- function(post_hoc_models_summaries) {
  .Deprecated("formatPHTResults")
  pkgs.require('purrr')
  
  
  # récupérer les différences entre les moyennes
  res_coeffs <- lapply(post_hoc_models_summaries, 
                       function(x) x$test$coefficients)
  logFC <- data.frame(purrr::map_dfr(res_coeffs, cbind),
                      row.names = names(post_hoc_models_summaries[[1]]$test$coefficients)
  )
  logFC <- as.data.frame(t(logFC))
  # extract raw p-values (non-adjusted)
  res_pvals <- lapply(post_hoc_models_summaries, function(x) x$test$pvalues)
  pvals <- data.frame(purrr::map_dfr(res_pvals, cbind),
                      row.names = names(post_hoc_models_summaries[[1]]$test$coefficients)
  )
  pvals <- as.data.frame(t(pvals))
  res <- list(
    "logFC" = logFC,
    "P_Value" = pvals
  )
  # formatting of column names for consistency with the limma and t-test code
  colnames(res$logFC) <- stringr::str_replace(
    colnames(res$logFC), " - ", "_vs_")
  colnames(res$P_Value) <- stringr::str_replace(
    colnames(res$P_Value), " - ", "_vs_")
  colnames(res$logFC) <- stringr::str_c(
    colnames(res$logFC), "_logFC")
  colnames(res$P_Value) <- stringr::str_c(
    colnames(res$P_Value), "_pval")
  
  return(res)
}





#' Post-hoc tests for classic 1-way ANOVA
#'
#' @description This function allows to compute a post-hoc test after a 1-way
#' ANOVA analysis. It expects as input an object obtained with the function
#' \code{classic1wayAnova}. The second parameter allows to choose between 2
#' different post-hoc tests: the Tukey Honest Significant Differences
#' (specified as "TukeyHSD") and the Dunnett test (specified as "Dunnett").
#'
#'
#' @param aov_fits a list containing aov fitted model objects
#'
#' @param post_hoc_test a character string indicating which post-hoc test to
#' use. Possible values are "TukeyHSD" or "Dunnett". See details for what to
#' choose according to your experimental design.
#'
#' @details
#' This is a function allowing to realise post-hoc tests for a set of
#' proteins/peptides for which a classic 1-way anova has been performed with
#' the function \code{classic1wayAnova}. Two types of tests are currently
#' available: The Tukey HSD's test and the Dunnett's test. Default is Tukey's
#' test.
#' The Tukey HSD's test compares all possible pairs of means, and is based on a
#' studentized range distribution. Here is used the \code{TukeyHSD()} function,
#' which can be applied to balanced designs (same number of samples in each
#' group), but also to midly unbalanced designs.
#' The Dunnett's test compares a single control group to all other groups.
#' Make sure the factor levels are properly ordered.
#'
#' @return a list of 2 dataframes: first one called "LogFC" contains
#' all pairwise comparisons logFC values (one column for one comparison) for
#' each analysed feature; The second one named "P_Value" contains the
#' corresponding pvalues.
#'
#' @author Hélène Borges
#'
#' @examples
#' \dontrun{examples/ex_postHocTest.R}
#'
#' @export
#'
#'
postHocTest <- function(aov_fits, post_hoc_test = "TukeyHSD") {
  .Deprecated("The other functions present in the file anova_analysis.R")
  pkgs.require('multcomp')
  
  
  if (post_hoc_test == "TukeyHSD") {
    # use of adjusted("none") to obtain raw p-values (and not adjusted ones)
    tukey_models_summaries <- lapply(
      aov_fits,
      function(x) {
        summary(multcomp::glht(x, 
                               linfct = multcomp::mcp(conditions = "Tukey")),
                test = multcomp::adjusted("none")
        )
      }
    )
    res <- formatPHResults(tukey_models_summaries)
  } else if (post_hoc_test == "Dunnett") {
    dunnett_models_summaries <- lapply(
      aov_fits,
      function(x) {
        summary(multcomp::glht(x, 
                               linfct = multcomp::mcp(conditions = "Dunnett")),
                test = multcomp::adjusted("none")
        )
      }
    )
    res <- formatPHResults(dunnett_models_summaries)
  } else {
    stop("Wrong post_hoc_test parameter. Please choose between TukeyHSD or
            Dunnett.")
  }
  return(res)
}
