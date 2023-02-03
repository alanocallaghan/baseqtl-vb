library("ggplot2")
library("latex2exp")
theme_set(theme_bw())

df <- readRDS("rds/noGT/constraints.rds")

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
ggsave("fig_1e-02/noGT/diag/constraints.pdf", width = 3, height = 4)
