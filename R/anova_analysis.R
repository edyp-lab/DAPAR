##' @title Function to perform a One-way Anova statistical test on a MsnBase dataset
##' @author Hélène Borges
##' @param currentLine The line currently treated from the quantitative data to
##'  perform the ANOVA
##' @param conditions The conditions represent the different classes of the
##' studied factor
##' @return A named vector containing all the different values of the aov model
classic1wayAnova <- function(current_line, conditions){
    # vector containing the protein/peptide intensities
    intensities <- unname(unlist(current_line))
    # anova sur ces deux vecteurs
    aov_test <- unlist(summary(aov(formula = intensities ~ conditions, data = NULL)))

    return(aov_test)

}

##' @title Wrapper for One-way Anova statistical test
##' @author Hélène Borges
##' @param obj An object of class \code{MSnSet}.
##' @param post_hoc a character string corresponding to the post-hoc test to
##' perform. Possible values are "TukeyHSD" and "Dunnett". See details of
##' \code{postHocTest()} function to choose the appropriate one.
##' @details This function allows to perform a 1-way Analysis of Variance. Also
##' computes the post-hoc tests if the \code{post_hoc} parameter is specified.
##' There are two possible post-hoc tests: the Tukey Honest Significant Differences
##' (specified as "TukeyHSD") and the Dunnett test (specified as "Dunnett").
##' @return A list of two dataframes. First one called "LogFC" contains
##' all pairwise comparisons logFC values (one column for one comparison) for
##' each analysed feature (Except in the case without post-hoc testing, for
##' which NAs are returned.); The second one named "P_Value" contains the
##' corresponding pvalues.
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' keepThat <- mvFilterGetIndices(obj, 'wholeMatrix', ncol(obj))
##' obj <- mvFilterFromIndices(obj, keepThat)
##' anovatest <- wrapperClassic1wayAnova(obj)
##' @seealso [postHocTest()]
wrapperClassic1wayAnova <- function(obj, post_hoc = NULL){
    qData <- Biobase::exprs(obj)
    sTab <- Biobase::pData(obj)
    if(post_hoc == "None" || is.null(post_hoc)){
        anova_tests <- as.data.frame(t(apply(qData,1,function(x) unlist(summary(classic1wayAnova(x,conditions=as.factor(sTab$Condition)))))))
        results <- dplyr::select(anova_tests, `Pr(>F)1`)
        to_return <- list("logFC" = data.frame("anova1way" = matrix(NA, nrow = nrow(results))),
                          "P_Value" = data.frame("anova1way" = p.adjust(results$`Pr(>F)1`, method = "fdr")))
    }else if(post_hoc == "TukeyHSD"){
        anova_tests <- t(apply(qData,1, classic1wayAnova, conditions=as.factor(sTab$Condition)))
        names(anova_tests) <- rownames(qData)
        to_return <- postHocTest(aov_fits = anova_tests, post_hoc_test = "TukeyHSD")

    }else if(post_hoc == "Dunnett"){
        anova_tests <- t(apply(qData,1, classic1wayAnova, conditions=as.factor(sTab$Condition)))
        names(anova_tests) <- rownames(qData)
        to_return <- postHocTest(aov_fits = anova_tests, post_hoc_test = "Dunnett")

    }else{
        stop("Wrong post_hoc parameter. Please choose between NULL, None,
             TukeyHSD or Dunnett.")
    }

    return(to_return)
}


##' Post-hoc tests for classic 1-way ANOVA
##'
##' @description This function allows to compute a post-hoc test after a 1-way
##' ANOVA analysis. It expects as input an object obtained with the function
##' \code{classic1wayAnova}. The second parameter allows to choose between 2
##' different post-hoc tests: the Tukey Honest Significant Differences
##' (specified as "TukeyHSD") and the Dunnett test (specified as "Dunnett").
##'
##' @param aov_object a list containing aov fitted model objects
##' @param post_hoc_test a character string indicating which post-hoc test to
##' use. Possible values are "TukeyHSD" or "Dunnett". See details for what to
##' choose according to your experimental design.
##' @details
##' This is a function allowing to realise post-hoc tests for a set of
##' proteins/peptides for which a classic 1-way anova has been performed with
##' the function \code{classic1wayAnova}. Two types of tests are currently
##' available: The Tukey HSD's test and the Dunnett's test. Default is Tukey's
##' test.
##' The Tukey HSD's test compares all possible pairs of means, and is based on a
##' studentized range distribution. Here is used the \code{TukeyHSD()} function,
##' which can be applied to balanced designs (same number of samples in each
##' group), but also to midly unbalanced designs.
##' The Dunnett's test compares a single control group to all other groups.
##' Make sure the factor levels are properly ordered.
##' @return a list of 2 dataframes: first one called "LogFC" contains
##' all pairwise comparisons logFC values (one column for one comparison) for
##' each analysed feature; The second one named "P_Value" contains the
##' corresponding pvalues.
##' @author Hélène Borges
postHocTest <- function(aov_fits, post_hoc_test = "TukeyHSD"){
    res <- NULL
    if(post_hoc_test == "TukeyHSD"){
        tukey_models <- lapply(aov_fits,TukeyHSD)
        # récupérer les différences entre les moyennes
        results_df_diffs <- lapply(tukey_models, function(x) data.frame(x$conditions,
                                                                        check.rows = TRUE, check.names = TRUE)$diff)
        logFC <- data.frame(purrr::map_dfr(results_df_diffs, cbind),
                            row.names = rownames(tukey_models[[1]]$conditions))
        logFC <- as.data.frame(t(logFC))
        # récupérer les pavalues ajustées
        results_df_pvals <- lapply(tukey_models, function(x) data.frame(x$conditions, check.rows = TRUE, check.names = TRUE)$p.adj)
        pvals <- data.frame(purrr::map_dfr(results_df_pvals, cbind),
                            row.names = rownames(tukey_models[[1]]$conditions))
        pvals <- as.data.frame(t(pvals))
        res <- list("logFC" = logFC,
                    "P_Value" = pvals)
    }else if(post_hoc_test == "Dunnett"){
        dunnett_models <- map(aov_fits,
                              multcomp::glht,
                              linfct = multcomp::mcp(conditions = "Dunnett"))
        res_summaries <- purrr::map(dunnett_models, summary)
        # on extrait les coefficients
        res_coefs <- purrr::map(res_summaries, function(x) x$test$coefficients)
        res_coefs <- data.frame(purrr::map_dfr(res_coefs, cbind),
                                row.names = names(res_summaries[[1]]$test$coefficients))
        res_coefs <- as.data.frame(t(res_coefs))
        # On extrait les pvalues
        res_pvalues <- purrr::map(res_summaries, function(x) x$test$pvalues)
        res_pvalues <- data.frame(purrr::map_dfr(res_pvalues, cbind),
                                  row.names = names(res_summaries[[1]]$test$coefficients))
        res_pvalues <- as.data.frame(t(res_pvalues))
        res <- list("logFC" = res_coefs,
                    "P_Value" = res_pvalues)
    }else{
        stop("Wrong post_hoc_test parameter. Please choose between TukeyHSD or
             Dunnett.")
    }
    return(res)
}
