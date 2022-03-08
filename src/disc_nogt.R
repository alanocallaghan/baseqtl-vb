library("baseqtl")
library("argparse")
library("rstan")
library("dplyr")

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--method",
    type = "character"
)
args <- parser$parse_args()

if (is.null(method <- args[["method"]])) {
    method <- "vb"
}
optimizing <- function(...) {
    rstan::optimizing(..., draws = 1000)
}
vb <- function(...) {
    rstan::vb(..., tol_rel_obj = 1e-3, iter = 20000)
}
sampling <- function(...) {
    rstan::sampling(..., chains = 1, open_progress = FALSE)
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

covariates <- list(
    normal_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/normal_skin_scaled_gc_libsize.rds"),
    Psoriasis_skin = readRDS("/home/abo27/rds/rds-mrc-bsu/ev250/alan/data/Psoriasis_skin_scaled_gc_libsize.rds")
)
disc <- readRDS("rds/noGT_discrepancies_vb.rds")


fit_stan <- function(i) {
    line <- disc[i, ]
    gene <- line[["gene"]]
    files <- c(
        normal_skin = sprintf("%s/refbias.%s.normal_skin.noGT.stan.input.rds", dir, gene),
        Psoriasis_skin = sprintf("%s/refbias.%s.Psoriasis_skin.noGT.stan.input.rds", dir, gene)
    )
    condition <- line[["condition"]]
    gene_data_cond <- readRDS(files[[condition]])
    snp <- line[["snp"]]
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
        draws <- 
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
disc_res <- lapply(
    1:nrow(disc),
    function(i) {
        fit_stan(i)
    }
)

disc1 <- disc %>% arrange(gene, snp, condition)

df <- do.call(rbind, disc_res)
df <- df %>% arrange(gene, snp, condition)

df$test <- paste(df$gene, df$snp, df$condition, sep = "_")
disc1$test <- paste(disc1$gene, disc1$snp, disc1$condition, sep = "_")


g <- ggplot() +
    aes(disc$vb, df$mean, colour = disc$discrepancy) +
    geom_point()
ggsave("tmp.png")


g <- ggplot() +
    aes(disc$hmc, df$mean, colour = disc$discrepancy) +
    geom_point()
ggsave("tmp1.png")
