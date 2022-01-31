library("baseqtl")
library("argparse")
library("rstan")

# probs <- c(0.005, 0.025, 0.25, 0.50, 0.75, 0.975, 0.995)
probs <- seq(0.005, 0.995, by = 0.005)
dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

files <- list.files(dir)
files <- grep("rbias", files, value=TRUE)
stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))

ls_mat <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")

gdd <- lapply(genes, function(gene) {
    gene_datas <- readRDS(
        sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
    )
    lapply(gene_datas, function(x) {
        out <- in.neg.beta.prob.eff2(x)
        out$cov <- cbind(out$cov, ls_mat[, gene])
        out$k <- 2
        out$aveP <- c(0, 0)
        out$sdP <- c(0.0309, 0.3479)
        out$mixP <- c(0.97359164, 0.02640836)
        out
    })
})
names(gdd) <- genes
saveRDS(gdd, "rds/gt_data.rds")


fit_stan <- function(i) {
    gene <- genes[[i]]
    cat(i, "/", length(genes), "\n")
    gene_datas <- readRDS(
        sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
    )
    snps_data <- lapply(
        gene_datas,
        function(x) {
            out <- in.neg.beta.prob.eff2(x)
            out$cov <- cbind(out$cov, ls_mat[, gene])
            out <- in.neg.beta.prob.eff2(x)
            out$k <- 2
            out$aveP <- c(0, 0)
            out$sdP <- c(0.0309, 0.3479)
            out$mixP <- c(0.97359164, 0.02640836)
            capture.output(
                time <- system.time(
                    post <- sampling(
                        baseqtl:::stanmodels$GT_nb, data = out,
                        chains = 1,
                        open_progress = FALSE
                    )
                )
            )
            tab <- rstan::summary(
                post,
                pars = "bj",
                probs = probs
            )$summary
            tab <- as.data.frame(tab)
            tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
            tab$gene <- gene
            tab$time <- time[["elapsed"]]
            tab
        }
    )
    names(snps_data) <- names(gene_datas)
    snps_data
}
all_stan_res <- parallel::mclapply(
    seq_along(genes),
    function(i) {
        file <- sprintf("rds/gt_nb_res/%s.rds", genes[[i]])
        if (file.exists(file)) {
            readRDS(file)
        } else {
            x <- tryCatch(
                fit_stan(i),
                error = function(e) list()
            )
            saveRDS(x, file)
        }
    }, mc.cores = 16
)

names(all_stan_res) <- genes

saveRDS(all_stan_res, "rds/gt_nb_res/all.rds")
