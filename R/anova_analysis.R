##' @title Function to perform a One-way Anova statistical test on a MsnBase dataset
##' @author Hélène Borges
##' @param currentLine The line currently treated from the quantitative data to
##'  perform the ANOVA
##' @param conditions The conditions represent the different classes of the
##' studied factor
##' @return A named vector containing al the differents values of the aov model
classic1wayAnova <- function(current_line, conditions){
    # vecteur contenant les abondances
    abondances <- unname(unlist(current_line))
    # anova sur ces deux vecteurs
    aov_test <- unlist(summary(aov(formula = abondances ~ conditions, data = NULL)))

    return(aov_test)

}

##' @title Wrapper for One-way Anova statistical test
##' @author Hélène Borges
##' @param obj An object of class \code{MSnSet}.
##' @return A dataframe of two columns. The first contains the p-values and the
##' second the adjusted p-values obtained withe Benjamini-Hochberg (FDR) method.
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept[1:1000]
##' keepThat <- mvFilterGetIndices(obj, 'wholeMatrix', ncol(obj))
##' obj <- mvFilterFromIndices(obj, keepThat)
##' anovatest <- wrapper.ow.aov(obj)
wrapperClassic1wayAnova <- function(obj){
    qData <- Biobase::exprs(obj)
    sTab <- Biobase::pData(obj)
    anova_test <- as.data.frame(t(apply(qData,1,classic1wayAnova,conditions=as.factor(sTab$Condition))))
    anova_test <- dplyr::select(anova_test, `Pr(>F)1`)
    anova_test$adjPVal <- p.adjust(anova_test$`Pr(>F)1`, method = "fdr")
    return(anova_test)
}
