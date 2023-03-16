library("ggplot2")
library("latex2exp")
library("argparse")
theme_set(theme_bw())

parser <- ArgumentParser()
parser$add_argument(
    "-m", "--model",
    default = "noGT",
    type = "character"
)
args <- parser$parse_args()

model <- args$model

df <- readRDS(sprintf("rds/%s/constraints.rds", model))

df <- df[df$model %in% c("master", "cons"), ]
df$model <- factor(df$model, levels = c("master", "cons"))
levels(df$model) <- c("No constraints", "With constraints")
df$method <- factor(df$method, levels = c("sampling", "vb"))
levels(df$method) <- c("HMC", "ADVI")
df$mean <- unlist(df$mean)
g <- ggplot(df, aes(x = method, y = mean)) +
    geom_jitter(alpha = 0.5) +
    facet_wrap(~model) +
    labs(x = "Inference method", y = TeX("Estimated $\\beta_{aFC}"))
ggsave(sprintf("fig_1e-02/%s/diag/constraints.pdf", model), width = 3, height = 2)
