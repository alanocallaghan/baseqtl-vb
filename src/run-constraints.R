library("baseqtl")
library("rstan")
library("argparse")
rstan_options(auto_write = TRUE)
source("src/functions.R")


parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "noGT",
    type = "character"
)
args <- parser$parse_args()

model <- args$model


if (model == "noGT") {
    files <- list.files(
        "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/",
        pattern = "refbias.*stan.input.rds",
        full.names = TRUE
    )
    datas <- lapply(files,
        function(file) {
            list(
                file = file,
                data = readRDS(file),
                gene = gsub(".*(ENSG\\d+)\\..*", "\\1", file),
                cond = gsub(".*(normal_skin|Psoriasis_skin).*", "\\1", file)
            )
        }
    )

    nulls <-  lapply(datas, function(x) {
        sapply(x$data, function(y) is.null(y$ase$ai0))
    })
    proportion_no_refbias <- sapply(nulls, sum) / sapply(nulls, length)
    ind_no_refbias <- which(proportion_no_refbias == 1)
    genes_no_refbias <- unique(sapply(datas, function(x) x$gene)[ind_no_refbias])
    gene_ind <- ind_no_refbias[[1]]
    gene <- datas[[gene_ind]]$gene
    # "ENSG00000057149", "SERPINB3"
    condition <- datas[[gene_ind]]$cond
    gene_data <- datas[[gene_ind]]$data

    length(gene_data)
    # [1] 13
    ## 13 SNPs should be enough to see an effect

    covariates <- get_covariates("noGT")

    ## take a couple of examples
    ## run stan-opt-d version and regular with vb and hmc
    ## remove constraints from stan-opt-d version and rerun

    master <- stan_model("src/stan/master/noGT_nb_ase.stan")
    opt_d_cons <- stan_model("src/stan/opt-d-cons/noGT_nb_ase.stan")
    # opt_d_nocons <- stan_model("src/stan/opt-d-nocons/noGT_nb_ase.stan")


    mods <- list(master=master, cons=opt_d_cons, nocons=opt_d_nocons)

    set.seed(42)
    snps <- names(gene_data)
    res <- parallel::mclapply(snps, function(snp) {
        snp_in <- gene_data[[snp]]
        data <- in.neg.beta.noGT.eff2(
            snp_in,
            covar = covariates[[condition]][names(snp_in$NB$counts), gene, drop = FALSE]
        )
        data$K <- 2
        data$k <- 2
        data$aveP <- c(0, 0)
        data$sdP <- c(0.0309, 0.3479)
        data$mixP <- c(0.97359164, 0.02640836)

        res <- lapply(mods,
            function(mod) {
                list(hmc = mean(extract(sampling(mod, data))$bj))
            }
        )
        sampling_res <- as.data.frame(do.call(rbind, res))
        colnames(sampling_res) <- "mean"
        sampling_res$fails <- 0
        sampling_res$model <- rownames(sampling_res)
        sampling_res$snp <- snp
        sampling_res$method <- "sampling"
        vb_dfs <- replicate(
            20,
            {
                res <- lapply(mods,
                    function(mod) {
                        fails <- 0
                        while (TRUE) {
                            f <- try(rstan::vb(mod, data, tol_rel_obj = 1e-2, iter = 50000, output_samples = 4000))
                            if (inherits(f, "try-error")) {
                                fails <- fails + 1
                            } else {
                                break
                            }
                        }
                        list(
                            advi = mean(extract(f)$bj),
                            fails = fails
                        )
                    }
                )
                vb_res <- as.data.frame(do.call(rbind, res))
                colnames(vb_res) <- c("mean", "fails")
                vb_res$model <- rownames(vb_res)
                vb_res$snp <- snp
                vb_res$method <- "vb"
                vb_res
            },
            simplify = FALSE
        )
        vb_df <- do.call(rbind, vb_dfs)
        rbind(sampling_res, vb_df)
    }, mc.cores = 8)

    df <- do.call(rbind, res)
    saveRDS(df, "rds/noGT/constraints.rds")
} else {

    files <- list.files(
        "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT",
        # pattern = "refbias.*stan.input.rds",
        pattern = "rbias.*.GT.stan1.input.rds",
        full.names = TRUE
    )

    datas <- lapply(files,
        function(file) {
            list(
                file = file,
                data = readRDS(file),
                gene = gsub(".*(ENSG\\d+)\\..*", "\\1", file)
            )
        }
    )

    gene <- datas[[2]]$gene
    # "ENSG00000025708", "TYMP"
    condition <- datas[[2]]$cond
    gene_data <- datas[[2]]$data

    length(gene_data)
    # [1] 33

    covariates <- get_covariates("GT")

    ## take a couple of examples
    ## run stan-opt-d version and regular with vb and hmc
    ## remove constraints from stan-opt-d version and rerun

    master <- stan_model("src/stan/master/GT_nb_ase.stan")
    opt_d_cons <- stan_model("src/stan/opt-d-cons/GT_nb_ase.stan")


    mods <- list(master=master, cons=opt_d_cons)

    set.seed(42)
    snps <- names(gene_data)
    res <- parallel::mclapply(snps, function(snp) {
        snp_in <- gene_data[[snp]]

        data <- in.neg.beta.prob.eff2(snp_in)
        data$cov <- cbind(data$cov, covariates[, gene])
        data$K <- 2
        data$k <- 2
        data$aveP <- c(0, 0)
        data$sdP <- c(0.0309, 0.3479)
        data$mixP <- c(0.97359164, 0.02640836)

        res <- lapply(mods,
            function(mod) {
                list(hmc = mean(extract(sampling(mod, data))$bj))
            }
        )
        sampling_res <- as.data.frame(do.call(rbind, res))
        colnames(sampling_res) <- "mean"
        sampling_res$fails <- 0
        sampling_res$model <- rownames(sampling_res)
        sampling_res$snp <- snp
        sampling_res$method <- "sampling"
        vb_dfs <- replicate(
            20,
            {
                res <- lapply(mods,
                    function(mod) {
                        fails <- 0
                        while (TRUE) {
                            f <- try(rstan::vb(mod, data, tol_rel_obj = 1e-2, iter = 50000, output_samples = 4000))
                            if (inherits(f, "try-error")) {
                                fails <- fails + 1
                            } else {
                                break
                            }
                        }
                        list(
                            advi = mean(extract(f)$bj),
                            fails = fails
                        )
                    }
                )
                vb_res <- as.data.frame(do.call(rbind, res))
                colnames(vb_res) <- c("mean", "fails")
                vb_res$model <- rownames(vb_res)
                vb_res$snp <- snp
                vb_res$method <- "vb"
                vb_res
            },
            simplify = FALSE
        )
        vb_df <- do.call(rbind, vb_dfs)
        rbind(sampling_res, vb_df)
    }, mc.cores = 8)

    df <- do.call(rbind, res)
    saveRDS(df, "rds/GT/constraints.rds")

}
