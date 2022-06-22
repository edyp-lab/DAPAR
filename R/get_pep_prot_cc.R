
#' @title Build the list of connex composant of the adjacency matrix
#'
#' @param X An adjacency matrix
#'
#' @return A list of CC
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept)
#' obj <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj, "Protein_group_IDs", FALSE)
#' ll <- get.pep.prot.cc(X)
#'
#' @export
#'
get.pep.prot.cc <- function(X) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Please install Matrix: BiocManager::install('Matrix')")
    }
    
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Please install igraph: BiocManager::install('igraph')")
    }
    
    if (!requireNamespace("graph", quietly = TRUE)) {
        stop("Please install graph: BiocManager::install('graph')")
    }
    
    if (is.null(X)) {
        warning("The adjacency matrix is empty")
        return()
    }

    # X <- as.matrix(X)
    p <- dim(X)[2] # Nb proteins
    q <- dim(X)[1] # Nb peptides


    multprot.cc <- singprot.cc <- multprot.cc.pep <- singprot.cc.pep <- NULL
    A <- B <- g <- NULL
    ### Adjacency matrix construction
    # boolean matrix product
    #print("Start computing boolean matrix product")
    A <- Matrix::crossprod(X, boolArith = TRUE)
    # A <- as.matrix(t(X) %*% X)
    # A <- as.matrix(t(X) %&% X)
    #print("End of computing boolean matrix product")
    # remove self-connecting edges
    diag(A) <- rep(0, p)
    # goes back to classical matrix format
    A <- matrix(as.numeric(A[, ]), ncol = p) 
    # reset pep and prot names
    colnames(A) <- rownames(A) <- colnames(X) 

    # proteins with no shared peptides
    # ie CC with only one protein
    SingleProt.CC.id <- which(rowSums(A) == 0)

    if (length(SingleProt.CC.id) > 0) {
        ### Peptides from single prot CCs
        singprot.cc <- as.list(names(SingleProt.CC.id))
        singprot.cc.pep <- list()
        for (i in seq_len(length(singprot.cc))) {
            peplist <- which(X[, singprot.cc[[i]]] != 0)
            singprot.cc.pep[[i]] <- names(peplist)
        }
    }


    # Test if there exist other CC than the single-prot ones
    otherCC.exists <- length(SingleProt.CC.id) < nrow(A)
    if (otherCC.exists) {
        # Remove all single-prot CC
        B <- A[-SingleProt.CC.id, -SingleProt.CC.id]

        ### Protein CCs
        # multprot.cc <- NULL
        # g <- graph::graphAM(B, edgemode='undirected', values=NA)
        # multprot.cc <- graph::connComp(as(g, 'graphNEL'))
        #
        # ### Peptides from multiple prot CCs
        # multprot.cc.pep <- list()
        # for(i in seq_len(length(multprot.cc))){
        #   protlist <- multprot.cc[[i]]
        #   subX <- as.matrix(X[,protlist])
        #   peplist <- which(rowSums(subX)!=0)
        #   multprot.cc.pep[[i]] <- names(peplist)
        # }

        multprot.cc <- NULL
        g2 <- igraph::graph.adjacency(B, mode = "undirected")
        cc.igraph <- igraph::components(g2)
        cc.id <- unique(cc.igraph$membership)
        multprot.cc <- lapply(
            cc.id,
            function(x) names(which(cc.igraph$membership == x))
        )

        multprot.cc.pep <- list()

        for (i in seq_len(length(multprot.cc))) {
            protlist <- multprot.cc[[i]]
            subX <- as.matrix(X[, protlist])
            peplist <- which(rowSums(subX) != 0)
            multprot.cc.pep[[i]] <- names(peplist)
        }
    }


    ### Merge results into a single list
    prot.cc <- c(multprot.cc, singprot.cc)
    pep.cc <- c(multprot.cc.pep, singprot.cc.pep)
    global.cc <- list()
    for (i in seq_len(length(prot.cc))) {
        prot <- prot.cc[[i]]
        pep <- pep.cc[[i]]
        tmp <- list(prot, pep)
        names(tmp) <- c("proteins", "peptides")
        global.cc[[i]] <- tmp
    }

    ### Clean memory and return result
    rm(
        A,
        B,
        g,
        multprot.cc,
        singprot.cc,
        multprot.cc.pep,
        singprot.cc.pep,
        prot.cc, pep.cc
    )
    gc()
    return(global.cc)
}





#' @title Jitter plot of CC
#'
#' @param list.of.cc List of cc such as returned by the function get.pep.prot.cc
#'
#' @return A plot
#'
#' @author Thomas Burger
#'
#' @examples
#' data(Exp1_R25_pept)
#' obj <- Exp1_R25_pept[seq_len(100)]
#' X <- BuildAdjacencyMatrix(obj, "Protein_group_IDs", TRUE)
#' ll <- get.pep.prot.cc(X)
#' plotJitter(ll)
#'
#' @export
#'
plotJitter <- function(list.of.cc = NULL) {
    if (is.null(list.of.cc)) {
        return()
    }

    #x <- length(list.of.cc) # number of CCs
    cc.summary <- sapply(list.of.cc, function(x) {
        c(length(x[[1]]), length(x[[2]]))
    })
    # cc.summary <- vapply(list.of.cc, 
    #     function(x) {c(length(x[[1]]), length(x[[2]]))},
    #     data.frame(2)
    #     )

    rownames(cc.summary) <- c("Nb_proteins", "Nb_peptides")
    colSums(cc.summary) # total amount of pep and prot in each CC
    colnames(cc.summary) <- seq_len(length(list.of.cc))
    cc.summary
    rowSums(cc.summary) # c(number of prot, number of pep)


    cc.summary <- as.data.frame(t(jitter(cc.summary)))
    plot(
        jitter(cc.summary[, 2]), 
        jitter(cc.summary[, 1]), 
        type = "p", 
        xlab = "#peptides in CC", 
        ylab = "#proteins in CC"
        )
}



#' @title Display a CC
#'
#' @param The.CC A cc (a list)
#'
#' @param X xxxxx
#'
#' @return A plot
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept)
#' obj <- Exp1_R25_pept[seq_len(100)]
#' X <- BuildAdjacencyMatrix(obj, "Protein_group_IDs", FALSE)
#' ll <- get.pep.prot.cc(X)
#' g <- buildGraph(ll[[1]], X)
#'
#' @export
#'
buildGraph <- function(The.CC, X) {
    subX <- as.matrix(X[The.CC$peptides, The.CC$proteins])
    nb.prot <- length(The.CC$proteins)
    nb.pep <- length(The.CC$peptides)
    nb.pep.shared <- length(which(rowSums(subX) > 1))
    nb.pep.spec <- length(which(rowSums(subX) == 1))
    nb.total <- nb.prot + nb.pep
    edge.list <- as.data.frame(which(subX == 1, arr.ind = TRUE))

    def.grp <- c(rep("shared.peptide", nb.pep), rep("protein", nb.prot))
    def.grp[which(rowSums(subX) == 1)] <- "spec.peptide"


    nodes <- data.frame(
        id = seq_len(nb.total),
        group = def.grp,
        label = c(rownames(subX), colnames(subX)),
        title = paste0("<p>", seq_len(nb.total), "<br>Tooltip !</p>"),
        size = c(rep(10, nb.pep), rep(20, nb.prot)),
        stringsAsFactors = FALSE
    )
    edges <- data.frame(
        from = c(edge.list$row),
        to = c(edge.list$col + nb.pep),
        stringsAsFactors = FALSE
    )

    return(
        list(
            nodes = nodes, 
            edges = edges
            )
        )
}


#' @title Display a CC
#'
#' @param g A cc (a list)
#'
#' @param layout xxxxx
#'
#' @param obj xxx
#'
#' @param prot.tooltip xxx
#'
#' @param pept.tooltip xxx
#'
#' @return A plot
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept)
#' obj <- Exp1_R25_pept[seq_len(100)]
#' X <- BuildAdjacencyMatrix(obj, "Protein_group_IDs", FALSE)
#' ll <- get.pep.prot.cc(X)
#' g <- buildGraph(ll[[1]], X)
#' display.CC.visNet(g)
#'
#' @export
#'
#'
#'
display.CC.visNet <- function(
    g,
    layout = layout_nicely,
    obj = NULL,
    prot.tooltip = NULL,
    pept.tooltip = NULL) {
    
    if (!requireNamespace("visNetwork", quietly = TRUE)) {
        stop("Please install visNetwork: BiocManager::install('visNetwork')")
    }

    col.prot <- "#ECB57C"
    col.spec <- "#5CA3F7"
    col.shared <- "#0EA513"


    visNetwork::visNetwork(g$nodes, g$edges, width = "100%", 
        height = "100%") %>%
        visNetwork::visNodes(shape = "dot") %>% # square for all nodes
        visNetwork::visGroups(groupname = "spec.peptide", 
            color = col.spec) %>% # darkblue for group "A"
        visNetwork::visGroups(groupname = "shared.peptide", 
            color = col.shared) %>% # darkblue for group "A"
        visNetwork::visGroups(groupname = "protein", 
            color = col.prot, shape = "dot") %>%
        visNetwork::visOptions(highlightNearest = FALSE) %>%
        # visLegend()
        # visPhysics(stabilization = FALSE)%>%
        visNetwork::visEdges(color = "#A9A9A9", width = 2)
    # %>%
    # visIgraphLayout(layout = "layout_with_fr")
}



#' @title Display a a jitter plot for CC
#'
#' @param df xxxx
#'
#' @param clickFunction xxxx
#'
#' @return A plot
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @export
#' 
#' @examples 
#' data(Exp1_R25_pept)
#' obj <- Exp1_R25_pept[seq_len(100)]
#' X <- BuildAdjacencyMatrix(obj, "Protein_group_IDs", TRUE)
#' ll <- get.pep.prot.cc(X)
#' plotJitter_rCharts(ll)
#'
plotJitter_rCharts <- function(df, clickFunction = NULL) {
    xtitle <- "TO DO"

    if (is.null(clickFunction)) {
        clickFunction <-
            JS("function(event) 
                {
                Shiny.onInputChange('eventPointClicked', 
                [this.index]+'_'+ [this.series.name]);
                }")
    }

    i_tooltip <- which(startsWith(colnames(df), "tooltip"))
    txt_tooltip <- NULL
    
    if (length(i_tooltip) == 0){
        stop("There is no tooltip in the object.")
    }
    
    for (i in i_tooltip) {
        txt_tooltip <- paste(txt_tooltip, "<b>", gsub("tooltip_", "",
            colnames(df)[i],
            fixed = TRUE
        ),
        " </b>: {point.", colnames(df)[i], "} <br> ",
        sep = ""
        )
    }

    h1 <- highchart() %>%
        hc_add_series(data = df, type = "scatter") %>%
        my_hc_chart(zoomType = "xy", chartType = "scatter") %>%
        hc_legend(enabled = FALSE) %>%
        hc_yAxis(title = list(text = "Nb of proteins ic CC")) %>%
        hc_xAxis(title = list(text = "Nb of peptides ic CC")) %>%
        hc_tooltip(headerFormat = "", pointFormat = txt_tooltip) %>%
        hc_plotOptions(series = list(
            animation = list(duration = 100),
            cursor = "pointer",
            point = list(events = list(
                click = clickFunction
            ))
        )) %>%
        my_hc_ExportMenu(filename = "plotCC")


    return(h1)
}
