f_g <- function(g, bj) {
    ifelse(
        abs(g) == 1,
        log1p(exp(bj)) - log(2),
        ifelse(
            g == 2,
            bj,
            ifelse(
                g == 0,
                0,
                NA
            )
        )
    )
}

gt_nb_gen <- function(
        N, cov = sim_cov(N), g = sim_g(N)
    ) {
    phi <- rgamma(1, shape = 1, rate = 0.1)
    alpha <- rnorm(1, 6, 4)
    beta <- rcauchy(1, 0, 2.5)
    which <- sample(1:4, 1)
    sizes <- c(0, 0.41, 0.55, 0.62)
    bj <- (sizes[which] + rnorm(1, sd = 0.01)) * sample(c(1, -1), 1)
    mu <- exp(f_g(g, bj) + alpha + (cov[, 3] * beta))
    Y <- rnbinom(N, mu = mu, size = phi)
    list(
        Y = Y,
        g = g,
        cov = cov,
        bj = bj,
        phi = phi,
        alpha = alpha,
        mu = mu,
        null = which == 1,
        beta = beta
    )
}

gt_nb_ase_gen <- function(
        N, cov = sim_cov(N), g = sim_g(N)
    ) {
    phi <- rgamma(1, shape = 1, rate = 0.1)
    alpha <- rnorm(1, 6, 4)
    beta <- rcauchy(1, 0, 2.5)
    which <- sample(1:4, 1)
    sizes <- c(0, 0.41, 0.55, 0.62)
    bj <- (sizes[which] + rnorm(1, sd = 0.01)) * sample(c(1, -1), 1)
    mu <- exp(f_g(g, bj) + alpha + (cov[, 3] * beta))
    Y <- rnbinom(N, mu = mu, size = phi)
    pi_i <- ifelse(g == 1, logistic(bj), 0.5)
    t <- Y / 10
    p_i <- rbeta(length(pi_i), shape1 = pi_i, shape2 = 5)
    s_i <- rbinom(length(pi_i), size = round(t), prob = p_i)
    list(
        Y = Y,
        g = g,
        cov = cov,
        bj = bj,
        phi = phi,
        alpha = alpha,
        mu = mu,
        m = round(t),
        p_i = p_i,
        n = s_i,
        null = which == 1,
        beta = beta
    )
}

logit <- function(x) {
    log(p / (1 - p))
}

logistic <- function(x) {
    1 / (1 + exp(-x))
}

sim_cov <- function(N) {
    m <- matrix(rnorm(N * 3), ncol = 3)
    m[, 1:2] <- 1
    m
}

sim_g <- function(N) {
    sample(0:2, N, replace = TRUE)
}
