library("baseqtl")
library("argparse")
library("rstan")

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--method",
    type = "character"
)
args <- parser$parse_args()

if (is.null(method <- args[["method"]])) {
    method <- "optimizing"
}
optimizing <- function(...) {
    rstan::optimizing(..., draws = 1000)
}
fun <- match.fun(method)

# probs <- c(0.005, 0.025, 0.25, 0.50, 0.75, 0.975, 0.995)
probs <- seq(0.005, 0.995, by = 0.005)

dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"
files <- list.files(dir, full.names = TRUE)

ns <- files[grep("refbias.ENSG\\d+\\.normal_skin.noGT.stan.input.rds", files)]
ps <- files[grep("refbias.ENSG\\d+\\.Psoriasis_skin.noGT.stan.input.rds", files)]

ngenes <- gsub(".*(ENSG\\d+).*", "\\1", ns)
pgenes <- gsub(".*(ENSG\\d+).*", "\\1", ps)

genes <- intersect(ngenes, pgenes)

ncovariates <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/normal_skin_scaled_gc_libsize.rds")
pcovariates <- readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/Psoriasis_skin_scaled_gc_libsize.rds")

gfiles <- lapply(
    genes,
    function(gene) {
        
    }
)

add <- c("N", "G", "A", "L") ## variables to add
same <- c("K", "k", "sdP", "aveP", "mixP") ## no change
bind <- "cov" ## to bind

vars <- c("ba", "bd", "bp", "bn")

fit_stan <- function(gene) {
    files <- c(
        normal = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, gene),
        psoriasis = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, gene)
    )
    no_all <- readRDS(files[[1]])
    pso_all <- readRDS(files[[2]])
    snps_both <- intersect(names(no_all), names(pso_all))
    lapply(snps_both,
        function(snp) {
            no <- no_all[[snp]]
            pso <- pso_all[[snp]]
            dn <- in.neg.beta.noGT.eff2(
                no,
                covar = ncovariates[names(no$NB$counts), gene, drop = FALSE]
            )
            dp <- in.neg.beta.noGT.eff2(
                pso,
                covar = pcovariates[names(pso$NB$counts), gene, drop = FALSE]
            )
            dp$k <- dn$k <- 3
            dp$aveP <- dn$aveP <- c(0, 0, 0)
            dp$sdP <- dn$sdP <- c(0.0436991990773286, 0.34926955206545, 0.4920048983496)
            dp$mixP <- dn$mixP <- c(-0.0460439385014068, -3.50655789731998, -4.19970507787993)

            out <- list()

            for (a in add) {
                out[[a]] <- dn[[a]] + dp[[a]]
            }
            for (s in same) {
                out[[s]] <- dn[[s]]
            }
            
            conc <- c("Y", "sNB", "gNB", "pNB", "gase", "m", "n", "pH", "ai0", "sdai0", "s", "h2g")
            if (is.null(dn$ai0) || is.null(dp$ai0)) {
                conc <- c("Y", "sNB", "gNB", "pNB", "gase", "m", "n", "pH", "s") ## to concatenate
            }
            for (x in conc) {
                out[[x]] <- c(dn[[x]], dp[[x]])
            }
            for (b in bind) {
                out[[b]] <- rbind(dn[[b]], dp[[b]])
            }

            l.sub <- list(dn, dp)
            seq <- c(0, sapply(l.sub, "[[", "A"))

            ## get cumsum, last element wont be used
            seq <- cumsum(seq)
            ## add number of ind with ASE from previous list to current list
            ## only to individuals with ASE counts (1 in col 1), add 0 to
            ## the first matrix.
            for (x in length(l.sub)) {
                l.sub[[x]][["ASEi"]][, 2][l.sub[[x]][["ASEi"]][, 1] == 1] <- l.sub[[x]][["ASEi"]][, 2][l.sub[[x]][["ASEi"]][, 1] == 1] + seq[x]
            }
            ## cbind and add to l
            out[["ASEi"]] <- Reduce(rbind, lapply(l.sub, "[[", "ASEi"))

            out[["I"]] <- rep(c(1, 0), sapply(l.sub, function(i) i[["N"]]))
            out[["IA"]] <- rep(c(1, 0), sapply(l.sub, function(i) i[["A"]]))

            if (is.null(dn$ai0) || is.null(dp$ai0)) {
                capture.output(
                    time <- system.time(
                        post <- fun(
                            baseqtl:::stanmodels$noGT2T_nb_ase_refbias, data = out
                        )
                    )
                )
            } else {
                capture.output(
                time <- system.time(
                    post <- fun(
                        baseqtl:::stanmodels$noGT2T_nb_ase, data = out
                    )
                )
            )
            }
            if (method == "optimizing") {
                tab <- posterior::summarise_draws(
                    post$theta_tilde,
                    mean,
                    sd,
                    ~quantile(.x, probs = probs)
                )
                tab <- tab[tab$variable %in% vars, ]
            } else {
                tab <- rstan::summary(
                    post,
                    pars = vars,
                    probs = probs
                )$summary
            }
            tab <- as.data.frame(tab)
            tab$null.99 <- sign(tab$"0.5%") == sign(tab$"99.5%")
            tab$gene <- gene
            tab$snp <- snp
            tab$time <- time[["elapsed"]]
            tab
        }
    )
}
all_stan_res <- parallel::mclapply(
    seq_along(genes),
    function(i) {
        cat(i, "/", length(genes), "\n")
        gene <- genes[[i]]
        file <- sprintf("rds/noGT/%s/%s.rds", method, gene)
        if (file.exists(file)) {
            readRDS(file)
        } else {
            x <- tryCatch(
                fit_stan(gene),
                error = function(e) list()
            )
            saveRDS(x, file)
        }
    }, mc.cores = 16
)
names(all_stan_res) <- genes

saveRDS(all_stan_res, sprintf("rds/noGT/%s/all.rds", method))
