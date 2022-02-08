library("ggplot2")
library("rjson")

pymc <- list.files("pymc_samples/")
pymc_df <- data.frame(
    gene = gsub("(.*)_(.*).json", "\\1", pymc),
    snp = gsub("(.*)_(.*).json", "\\2", pymc),
    test = gsub("(.*_.*).json", "\\1", pymc)
)
pymc_df$file <- list.files("pymc_samples/", full.names = TRUE)

pymc_df$mean <- NA
pymc_df$time <- NA
pymc_df$sd <- NA
for (i in seq_along(pymc_df$file)) {
    cat(i, "/", nrow(pymc_df), "\n")
    d <- fromJSON(file = pymc_df$file[[i]])
    pymc_df$mean[[i]] <- d$samples$location
    pymc_df$sd[[i]] <- d$samples$scale
    pymc_df$time[[i]] <- d$time
}



pyro <- fromJSON(file = "pyro_params.json")

pyro_df <- data.frame(
    gene = gsub("(.*)_(.*)", "\\1", names(pyro)),
    snp = gsub("(.*)_(.*)", "\\2", names(pyro)),
    test = names(pyro)
)

pyro_df$mean <- NA
pyro_df$time <- NA
pyro_df$sd <- NA
for (i in seq_len(nrow(pyro_df))) {
    cat(i, "/", nrow(pymc_df), "\n")
    pyro_df$mean[[i]] <- pyro[[i]]$params$location
    pyro_df$sd[[i]] <- pyro[[i]]$params$scale
    pyro_df$time[[i]] <- pyro[[i]]$time
}


stan_files <- list.files("rds/gt_nb_res/", pattern = "ENSG*", full.names = TRUE)
genes <- gsub(".*(ENSG.*).rds", "\\1", stan_files)
stan_list <- lapply(seq_along(stan_files),
    function(i) {
        f <- stan_files[[i]]
        x <- readRDS(f)
        df <- do.call(rbind, x)
        df$snp <- names(x)
        df$gene <- genes[[i]]
        df
    }
)
stan_df <- do.call(rbind, stan_list)
stan_df$test <- paste(stan_df$gene, stan_df$snp, sep = "_")


pyro_pymc <- merge(pymc_df, pyro_df, by = "test", suffix = c(".pymc", ".pyro"))
stan_pyro_pymc <- merge(stan_df, pyro_pymc, by = "test", suffix = c("stan", ".pyro"))


g <- ggplot(stan_pyro_pymc) +
    aes(mean, mean.pyro) +
    geom_point()
ggsave("stan-pyro.png")

g <- ggplot(stan_pyro_pymc) +
    aes(mean, mean.pymc) +
    geom_point()
ggsave("stan-pymc.png")

g <- ggplot(stan_pyro_pymc) +
    aes(mean.pymc, mean.pyro) +
    geom_point()
ggsave("pymc-pyro.png")

g <- ggplot(pyro_pymc) +
    geom_density(aes(x = time.pyro, colour = "pyro")) +
    geom_density(aes(x = time.pymc, colour = "pymc")) +
    scale_x_log10()
ggsave("time-pymc-pyro.png")
