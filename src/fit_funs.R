if (!exists("method")) {
    method <- "vb"
}
fun <- match.fun(method)

probs <- seq(0.005, 0.995, by = 0.005)
if (model == "GT") {
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

    files <- list.files(dir)
    files <- grep("rbias", files, value=TRUE)
    stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
    genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))

    ls_mat <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")

    fit_stan <- function(gene, snp, init = "random") {
        gene_datas <- readRDS(
            sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
        )
        x <- gene_datas[[snp]]
        data <- in.neg.beta.prob.eff2(x)
        data$cov <- cbind(data$cov, ls_mat[, gene])
        data$K <- 2
        data$k <- 2
        data$aveP <- c(0, 0)
        data$sdP <- c(0.0309, 0.3479)
        data$mixP <- c(0.97359164, 0.02640836)
        if (init == "opt") {
            opt_pars <- optimizing(
                baseqtl:::stanmodels$GT_nb_ase_refbias, data = data
            )
            init <- extract_params(opt_pars)
            if (method == "sampling") {
                init <- rep(init, 4)
            }
        } else if (init == "fixed") {
            init <- list(
                betas = c(5, 0),
                bj = 0,
                phi = 10,
                theta = 10,
                rai0 = rep(0, data$L)
            )
        }
        while (TRUE) {
            t0 <- proc.time()
            capture.output(
                post <- try(
                    fun(
                        baseqtl:::stanmodels$GT_nb_ase_refbias, data = data,
                        init = init
                    )
                )
            )
            time <- proc.time() - t0
            if (!inherits(post, "try-error")) break
        }
        if (method == "optimizing") {
            tab <- posterior::summarise_draws(
                post$theta_tilde,
                mean,
                sd,
                ~quantile(.x, probs = probs)
            )
            tab <- tab[tab$variable == "bj", ]
        } else {
            tab <- rstan::summary(
                post,
                pars = "bj",
                probs = probs
            )$summary
        }
        tab <- as.data.frame(tab)
        tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
        tab$gene <- gene
        tab$time <- time[["elapsed"]]
        tab$snp <- snp
        tab
    }
} else {
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"
    files <- list.files(dir, full.names = TRUE)

    ns <- files[grep("refbias.ENSG\\d+\\.normal_skin.noGT.stan.input.rds", files)]
    ps <- files[grep("refbias.ENSG\\d+\\.Psoriasis_skin.noGT.stan.input.rds", files)]

    ngenes <- gsub(".*(ENSG\\d+).*", "\\1", ns)
    pgenes <- gsub(".*(ENSG\\d+).*", "\\1", ps)

    genes <- intersect(ngenes, pgenes)

    covariates <- list(
        normal_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/normal_skin_scaled_gc_libsize.rds"),
        Psoriasis_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/Psoriasis_skin_scaled_gc_libsize.rds")
    )

    fit_stan <- function(gene, snp, init = "random") {
        files <- c(
            normal_skin = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, gene),
            Psoriasis_skin = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, gene)
        )
        out <- lapply(names(files),
            function(condition) {
                gene_data_cond <- readRDS(files[[condition]])
                snp_in <- gene_data_cond[[snp]]
                data <- in.neg.beta.noGT.eff2(
                    snp_in,
                    covar = covariates[[condition]][names(snp_in$NB$counts), gene, drop = FALSE]
                )
                data$k <- 3
                data$aveP <- c(0, 0, 0)
                data$sdP <- c(0.0436991990773286, 0.34926955206545, 0.4920048983496)
                data$mixP <- c(-0.0460439385014068, -3.50655789731998, -4.19970507787993)

                model <- baseqtl:::stanmodels$noGT_nb_ase
                if (!is.null(data$ai0)) {
                    model <- baseqtl:::stanmodels$noGT_nb_ase_refbias
                }

                if (init == "opt") {
                    opt_pars <- optimizing(model, data = data)
                    init <- extract_params(opt_pars)
                    if (method == "sampling") {
                        init <- rep(init, 4)
                    }
                } else if (init == "fixed") {
                    init <- list(
                        betas = c(5, 0),
                        bj = 0,
                        phi = 10,
                        theta = 10,
                        rai0 = rep(0, data$L)
                    )
                }
                while (TRUE) {
                    t0 <- proc.time()
                    capture.output(
                        post <- try(
                            fun(
                                model, data = data,
                                init = init
                            )
                        )
                    )
                    time <- proc.time() - t0
                    if (!inherits(post, "try-error")) break
                }
                if (method == "optimizing") {
                    tab <- posterior::summarise_draws(
                        post$theta_tilde,
                        mean,
                        sd,
                        ~quantile(.x, probs = probs)
                    )
                    tab <- tab[tab$variable %in% "bj", ]
                } else {
                    tab <- rstan::summary(
                        post,
                        pars = "bj",
                        probs = probs
                    )$summary
                }
                tab <- as.data.frame(tab)
                tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
                tab$gene <- gene
                tab$snp <- snp
                tab$time <- time[["elapsed"]]
                tab$condition <- condition
                tab
            }
        )
        do.call(rbind, out)
    }
}


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
