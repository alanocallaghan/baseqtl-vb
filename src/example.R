library("baseqtl")
library("rstan")

gene <- "ENSG00000015475"
snp <- "17718119:G:T"

ls_mat <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")
dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

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
post <- rstan::sampling(
    baseqtl:::stanmodels$GT_nb_ase_refbias, data = out,
    chains = 4, open_progress = FALSE
)
tab <- rstan::summary(
    post,
    pars = "bj"
)$summary
tab <- as.data.frame(tab)
tab$gene <- gene
tab$snp <- snp



dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"
gene <- "ENSG00000002330"

covariates <- list(
    normal_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/normal_skin_scaled_gc_libsize.rds"),
    Psoriasis_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/Psoriasis_skin_scaled_gc_libsize.rds")
)

files <- c(
    normal_skin = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, gene),
    Psoriasis_skin = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, gene)
)

results <- lapply(
    names(files),
    function(condition) {
        gene_data_cond <- readRDS(files[[condition]])
        snp <- "63853645:C:T"
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
        post <- rstan::sampling(
            model,
            data = data,
            chains = 4,
            open_progress = FALSE
        )
        tab <- rstan::summary(
            post,
            pars = "bj"
        )$summary
        tab <- as.data.frame(tab)
        tab$gene <- gene
        tab$snp <- snp
        tab$condition <- condition
        tab
    }
)
