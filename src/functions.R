#' Plot the ELBO over iterations of the Stan variational Bayes algorithm ADVI.
#' 
#' @param x An object of class "scive".
#' @return A ggplot.
#' @examples
#' \dontrun{
#'  sce <- mockSCE()
#'  n <- scive(sce, method = "VB", tol_rel_obj = 0.1)
#'  plot_elbo(n)
#' }
#' @export
plot_elbo <- function(elbo) {
    if (is.null(elbo)) {
        stop(
            'Object does not contain ELBO information.\n',
            'They must have been created with `fun="vb"`.')
    }
    if (is.data.frame(elbo)) {
        ggplot(elbo, aes(x = iter, y = log(abs(ELBO)))) +
            geom_line()
    }
}

parse_elbo <- function(elbo) {
    elbo <- gsub("Chain 1:\\s+", "", elbo)
    normal <- grep("Drawing", elbo)
    abnormal <- grep("Informational", elbo)
    if (!length(normal)) {
        ## stop or warning????
        stop("Failed to parse ELBO")
    }
    ind_end <- normal - 2
    if (length(abnormal)) {
        ind_end <- abnormal - 1
    }
    elbo <- elbo[(grep("Begin stochastic", elbo) + 1):ind_end]
    ## This ain't quite right
    elbo <- gsub("MAY BE DIVERGING... INSPECT ELBO", "", elbo, fixed = TRUE)
    elbo[-c(1, grep("CONVERGED", elbo))] <- paste(
        elbo[-c(1, grep("CONVERGED", elbo))],
        "NOTCONVERGED"
    )
    ## If both mean and median elbo converge this ends up with CONVERGED CONVERGED
    ## and therefore another column
    elbo <- gsub("(MEDIAN |MEAN )?ELBO CONVERGED", "CONVERGED", elbo)
    elbo <- gsub("(\\s+CONVERGED){2}", " CONVERGED", elbo)
    elbo <- strsplit(elbo, "\\s+")
    elbo <- do.call(rbind, elbo)
    colnames(elbo) <- elbo[1, ]
    elbo <- elbo[-1, ]
    elbo <- as.data.frame(elbo, stringsAsFactors = FALSE)
    elbo[, 1:4] <- lapply(elbo[, 1:4], as.numeric)
    elbo
}
