if (!exists("method")) {
    method <- "vb"
}
fun <- match.fun(method)

if (model == "GT") {
    probs <- seq(0.005, 0.995, by = 0.005)
    dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

    files <- list.files(dir)
    files <- grep("rbias", files, value=TRUE)
    stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
    genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))

    ls_mat <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")

    fit_stan <- function(gene, snp) {
        gene_datas <- readRDS(
            sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
        )
        x <- gene_datas[[snp]]
        out <- in.neg.beta.prob.eff2(x)
        out$cov <- cbind(out$cov, ls_mat[, gene])
        out$K <- 2
        out$k <- 2
        out$aveP <- c(0, 0)
        out$sdP <- c(0.0309, 0.3479)
        out$mixP <- c(0.97359164, 0.02640836)
        capture.output(
            time <- system.time(
                post <- fun(
                    baseqtl:::stanmodels$GT_nb_ase_refbias, data = out
                )
            )
        )
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

    fit_stan <- function(gene, snp) {
        files <- c(
            normal_skin = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, gene),
            Psoriasis_skin = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, gene)
        )
        lapply(names(files),
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
                capture.output(
                    time <- system.time(
                        post <- fun(
                            model, data = data
                        )
                    )
                )
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
    }
}
