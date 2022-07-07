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
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(1000)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' anova_tests <- t(apply(Biobase::exprs(obj$new), 1, classic1wayAnova,
#'     conditions = as.factor(Biobase::pData(obj$new)$Condition)
#' ))
#'
#' @export
#'
classic1wayAnova <- function(current_line, conditions) {
    
    if (!requireNamespace("stats", quietly = TRUE)) {
        stop("Please install stats: BiocManager::install('stats')")
    }
    
    
    # vector containing the protein/peptide intensities
    intensities <- unname(unlist(current_line))
    # anova sur ces deux vecteurs
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
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(1000)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' anovatest <- wrapperClassic1wayAnova(obj$new)
#'
#' @seealso [postHocTest()]
#'
#' @export
#'
wrapperClassic1wayAnova <- function(
    obj, 
    with_post_hoc = "No", 
    post_hoc_test = "No") {
    
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Please install dplyr: BiocManager::install('dplyr')")
    }
    
    
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
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(1000)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' anova_tests <- t(apply(Biobase::exprs(obj$new), 1, classic1wayAnova,
#'     conditions = as.factor(Biobase::pData(obj$new)$Condition)
#' ))
#' names(anova_tests) <- rownames(Biobase::exprs(obj$new))
#' tms <- lapply(
#'     anova_tests,
#'     function(x) {
#'         summary(multcomp::glht(x,
#'             linfct = multcomp::mcp(conditions = "Tukey")
#'         ),
#'         test = multcomp::adjusted("none")
#'         )
#'     }
#' )
#' res <- formatPHResults(tms)
#'
#' @export
#'
formatPHResults <- function(post_hoc_models_summaries) {
    if (!requireNamespace("purrr", quietly = TRUE)) {
        stop("Please install purrr: BiocManager::install('purrr')")
    }

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
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(1000)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' anova_tests <- t(apply(Biobase::exprs(obj$new), 1, classic1wayAnova,
#'     conditions = as.factor(Biobase::pData(obj$new)$Condition)
#' ))
#' names(anova_tests) <- rownames(Biobase::exprs(obj$new))
#' pht <- postHocTest(aov_fits = anova_tests)
#'
#' @export
#'
#'
postHocTest <- function(aov_fits, post_hoc_test = "TukeyHSD") {
    if (!requireNamespace("multcomp", quietly = TRUE)) {
        stop("Please install multcomp: BiocManager::install('multcomp')")
    }

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
