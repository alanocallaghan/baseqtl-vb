library("baseqtl")
library("BiocParallel")
library("argparse")

probs <- seq(0.005, 0.995, by = 0.005)
dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"

files <- list.files(dir)
files <- grep("rbias", files, value = TRUE)
stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))

ls_mat <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/scaled_gc_libsize.rds")
register(MulticoreParam(workers = 8))

lms <- bplapply(seq_along(genes),
    function(i) {
        cat(i, "/", length(genes), "genes, \n")
        gene <- genes[[i]]
        gene_datas <- readRDS(
            sprintf("%s/rbias.%s.GT.stan1.input.rds", dir, gene)
        )
        mods <- lapply(
            seq_along(gene_datas),
            function(j) {
                # cat(j, "/", length(gene_datas), "snps \n")
                snp <- names(gene_datas)[[j]]
                x <- gene_datas[[snp]]
                out <- in.neg.beta.prob.eff2(x)
                out$cov <- cbind(out$cov, ls_mat[, gene])
                t0 <- proc.time()
                m1 <- cbind(
                    as.data.frame(
                        summary(lm(log(out$Y) ~ abs(out$g) + out$cov[, -(1:2)]))$coef
                    )[2, ],
                    method = "lm",
                    gene = gene,
                    snp = snp
                )
                t1 <- proc.time()
                m2 <- cbind(
                    as.data.frame(
                        summary(
                            MASS::glm.nb(out$Y ~ abs(out$g) + out$cov[, -(1:2)], link = "log")
                        )$coef
                    )[2, ],
                    method = "glm",
                    gene = gene,
                    snp = snp
                )
                t2 <- proc.time()
                colnames(m1) <- colnames(m2) <- c("coef", "stderr", "stat", "pval", "method", "gene", "snp")
                m1$time <- (t1 - t0)[["elapsed"]]
                m2$time <- (t2 - t1)[["elapsed"]]
                out <- rbind(m1, m2)
                out
            }
        )
        do.call(rbind, mods)
    },
    mc.cores = 8
)
df <- do.call(rbind, lms)
df$padj <- NA
df$padj[df$method == "lm"] <- p.adjust(df$pval[df$method == "lm"], method = "BH")
df$padj[df$method == "glm"] <- p.adjust(df$pval[df$method == "glm"], method = "BH")

saveRDS(df, "rds/GT/lm-filtering.rds")
