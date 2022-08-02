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

## this is the equivalent of rstan::extract but for `optimizing` results.
extract_params <- function(fit) {
    pars <- gsub("\\[(\\d+,)*\\d+\\]", "", names(fit$par))
    dims <- gsub("^[a-z0_]+\\[(.*)\\]", "\\1", names(fit$par))
    out <- lapply(
        unique(pars),
        function(p) {
            ind <- pars == p
            if (!all(grepl("[", names(fit$par)[ind], fixed = TRUE))) {
                return(fit$par[ind])
            }
            if (!all(grepl(",", names(fit$par)[ind], fixed = TRUE))) {
                return(fit$par[ind])
            }
            l <- strsplit(dims[ind], ",")
            d <- do.call(rbind, l)
            d <- matrix(as.numeric(d), ncol = ncol(d))
            array(
                fit$par[ind],
                dim = c(1, apply(d, 2, max))
            )
        }
    )
    setNames(out, unique(pars))
}

## some wrappers to define sensible defaults for stan
optimizing <- function(...) {
    while (TRUE) {
        f <- try(rstan::optimizing(...), silent = TRUE)
        if (!inherits(f, "try-error")) {
            return (f)
        }
    }
}

vb <- function(..., tol_rel_obj = 1e-2, n_iterations = 50000) {
    # rstan::vb(..., grad_samples = 5, elbo_samples = 1000, tol_rel_obj = 1e-3)
    while (TRUE) {
        f <- try(rstan::vb(..., tol_rel_obj = tol_rel_obj, iter = n_iterations), silent = TRUE)
        if (!inherits(f, "try-error")) {
            return (f)
        }
    }
}

sampling <- function(..., chains = 4) {
    rstan::sampling(..., chains = chains, open_progress = FALSE)
}

fit_stan_GT <- function(
        gene_data,
        gene,
        snp,
        covariates = get_covariates("GT"),
        init = "random",
        method = c("vb", "sampling", "optimizing"),
        tol = 1e-2,
        vars = "bj",
        model = baseqtl:::stanmodels$GT_nb_ase_refbias,
        probs = seq(0.005, 0.995, by = 0.005),
        summarise_posterior = TRUE,
        ...
    ) {
    method <- match.arg(method)
    fun <- match.fun(method)
    
    data <- in.neg.beta.prob.eff2(gene_data[[snp]])
    data$cov <- cbind(data$cov, covariates[, gene])
    data$K <- 2
    data$k <- 2
    data$aveP <- c(0, 0)
    data$sdP <- c(0.0309, 0.3479)
    data$mixP <- c(0.97359164, 0.02640836)

    if (identical(init, "opt")) {
        opt_pars <- optimizing(model, data = data)
        init <- extract_params(opt_pars)
    } else if (identical(init, "fixed")) {
        init <- list(
            betas = c(5, 0),
            bj = 0,
            phi = 10,
            theta = 10,
            rai0 = rep(0, data$L)
        )
    }
    if (is.list(init) && method == "sampling") {
        init <- replicate(4, init, simplify = FALSE)
    }
    t0 <- proc.time()
    capture.output(
        post <- fun(
            model,
            data = data,
            init = init,
            ...
        )
    )
    time <- proc.time() - t0
    if (!summarise_posterior) {
        return(post)
    }
    
    tab <- summarise_post(post, method, vars, probs)
    tab$gene <- gene
    tab$time <- time[["elapsed"]]
    tab$snp <- snp
    tab
}

summarise_post <- function(post, method, vars, probs) {
    if (method == "optimizing") {
        tab <- posterior::summarise_draws(
            post$theta_tilde,
            mean,
            sd,
            ~quantile(.x, probs = probs)
        )
        if (!is.null(vars)) {
            tab <- tab[tab$variable %in% vars, ]
        }
    } else {
        tab <- rstan::summary(
            post,
            probs = probs
        )$summary
        if (!is.null(vars)) {
            tab <- tab[rownames(tab) %in% vars, ]
        }
    }
    tab <- as.data.frame(tab)
    tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
    tab
}

fit_stan_noGT <- function(
        gene_data,
        gene,
        snp,
        covariates = get_covariates("noGT"),
        init = "random",
        method = c("vb", "sampling", "optimizing"),
        tol = 1e-2,
        vars = "bj",
        model = baseqtl:::stanmodels$noGT_nb_ase,
        probs = seq(0.005, 0.995, by = 0.005),
        summarise_posterior = TRUE,
        ...
    ) {

    method <- match.arg(method)
    fun <- match.fun(method)

    out <- lapply(names(gene_data),
        function(condition) {
            snp_in <- gene_data[[condition]][[snp]]
            data <- in.neg.beta.noGT.eff2(
                snp_in,
                covar = covariates[[condition]][names(snp_in$NB$counts), gene, drop = FALSE]
            )
            data$k <- 3
            data$aveP <- c(0, 0, 0)
            data$sdP <- c(0.0436991990773286, 0.34926955206545, 0.4920048983496)
            data$mixP <- c(-0.0460439385014068, -3.50655789731998, -4.19970507787993)

            ## ignoring this for the time being
            # model <- baseqtl:::stanmodels$noGT_nb_ase
            # if (!is.null(data$ai0)) {
            #     model <- baseqtl:::stanmodels$noGT_nb_ase_refbias
            # }

            if (identical(init, "opt")) {
                opt_pars <- optimizing(model, data = data, hessian = TRUE)
                init <- extract_params(opt_pars)
            } else if (identical(init, "fixed")) {
                init <- list(
                    betas = c(5, 0),
                    bj = 0,
                    phi = 10,
                    theta = 10,
                    rai0 = rep(0, data$L)
                )
            }
            if (is.list(init) && method == "sampling") {
                init <- replicate(4, init, simplify = FALSE)
            }
            t0 <- proc.time()
            capture.output(
                post <- fun(
                    model,
                    data = data,
                    init = init,
                    ...
                )
            )
            time <- proc.time() - t0
            if (!summarise_posterior) {
                return(post)
            }
            tab <- summarise_post(post, method, vars, probs)
            tab$gene <- gene
            tab$time <- time[["elapsed"]]
            tab$snp <- snp
            tab$condition <- condition
            tab
        }
    )
    names(out) <- names(gene_data)
    if (!summarise_posterior) {
        return(out)
    }
    do.call(rbind, out)
}


get_gene_data <- function(gene, model) {
    if (model == "GT") {
        dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"
        gene_data <- readRDS(
            sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
        )
    } else {
        dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"

        files <- c(
            normal_skin = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, gene),
            Psoriasis_skin = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, gene)
        )
        gene_data <- lapply(files, readRDS)
    }
    gene_data
}

get_snps <- function(gene_data) {
    if (identical(names(gene_data), c("normal_skin", "Psoriasis_skin"))) {
        snps_each <- lapply(gene_data, names)
        do.call(intersect, unname(snps_each))
    } else {
        names(gene_data)
    }
}


get_covariates <- function(model) {
    if (model == "GT") {
        readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")
    } else {
        list(
            normal_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/normal_skin_scaled_gc_libsize.rds"),
            Psoriasis_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/Psoriasis_skin_scaled_gc_libsize.rds")
        )
    }
}

fit_mod <- function(model, gene, snp, condition, method = "sampling", stanmodel = NULL, ...) {
    fun <- match.fun(method)
    if (model == "noGT") {
        gene_data <- get_gene_data(gene, model)
        covariates <- get_covariates("noGT")

        snp_in <- gene_data[[condition]][[snp]]
        data <- in.neg.beta.noGT.eff2(
            snp_in,
            covar = covariates[[condition]][names(snp_in$NB$counts), gene, drop = FALSE]
        )
        data$k <- 3
        data$aveP <- c(0, 0, 0)
        data$sdP <- c(0.0436991990773286, 0.34926955206545, 0.4920048983496)
        data$mixP <- c(-0.0460439385014068, -3.50655789731998, -4.19970507787993)

        stanmodel <- stanmodel %||% baseqtl:::stanmodels$noGT_nb_ase
        fit <- fun(
            stanmodel,
            data = data,
            ...
        )
    } else {
        gene_data <- get_gene_data(gene, model)
        covariates <- get_covariates("GT")

        data <- in.neg.beta.prob.eff2(gene_data[[snp]])
        data$cov <- cbind(data$cov, covariates[, gene])
        data$K <- 2
        data$k <- 2
        data$aveP <- c(0, 0)
        data$sdP <- c(0.0309, 0.3479)
        data$mixP <- c(0.97359164, 0.02640836)
        stanmodel <- stanmodel %||% baseqtl:::stanmodels$GT_nb_ase_refbias
        fit <- fun(
            stanmodel,
            data = data,
            ...
        )
    }
    fit
}

plot_with_legend_below <- function(
        ...,
        rel_heights = if (common_x) c(0.85, 0.05, 0.15) else c(0.85, 0.15),
        common_x = FALSE,
        nrow = 1,
        labels = "AUTO",
        align = "h",
        xlab
    ) {
    legends <- lapply(list(...), cowplot::get_legend)
    plots <- list(...)
    # if (!all(sapply(legends, function(x) identical(x$grobs, legends[[1]]$grobs)))) {
    #   stop("Different legends shouldn't be merged!")
    # }
    t <- theme(legend.position = "none")
    if (common_x) {
        t <- t + theme(axis.title.x = element_blank())
    }
    together <- cowplot::plot_grid(
        plotlist = lapply(plots, function(x) x + t),
        labels = labels,
        nrow = nrow,
        align = align
    )
    if (common_x) {
        if (missing(xlab)) {
            xlab <- rlang::as_label(plots[[1]]$mapping$x)
        }
        xlab <- grid::textGrob(
            xlab,
            gp = grid::gpar(fontsize = 12)
        )
        cowplot::plot_grid(
            together, xlab, legends[[1]], nrow = 3, rel_heights = rel_heights
        )
    } else {
        cowplot::plot_grid(
            together, legends[[1]], nrow = 2, rel_heights = rel_heights
        )
    }
}


# make a latex crosstab with suitable lines and labels
make_crosstab <- function(x, y, xn = "x", yn = "y", file = "", prop = FALSE, ...) {
    x <- factor(x)
    y <- factor(y)
    tab <- table(x, y)
    if (prop) {
        tab <- tab / sum(tab)
        tab <- round(tab, digits = 3)
    }
    tab <- as.data.frame(matrix(tab, ncol = nlevels(y)))
    colnames(tab) <- levels(y)
    tab <- cbind(" " = c(xn, rep(" ", length.out = nlevels(x) - 1)), " " = levels(x), tab)
    tab <- rbind(
        c("", "", yn, rep("", nlevels(y) - 1)),
        c(" ", " ", levels(y)),
        tab
    )
    xt <- xtable::xtable(
        tab,
        format = "latex",
        align = paste0("rrr|", paste(rep("r", nlevels(y)), collapse = "")),
        col.names = NULL,
        ...
    )
    print(xt,
        include.rownames = FALSE,
        include.colnames = FALSE,
        hline.after = 2,
        file = file
    )

}

expanded_range <- function(x, frac = 0.01, maxabs = NULL) {
    r <- range(x)
    d <- abs(r[[1]] - r[[2]])
    r[[1]] <- r[[1]] - (d * (1 + frac))
    r[[2]] <- r[[2]] + (d * (1 + frac))
    if (!is.null(maxabs)) {
        r[abs(r) > maxabs] <- maxabs * sign(r)
    }
    r
}

`%||%` <- function(a, b) {
    if (is.null(a)) b else a
}
