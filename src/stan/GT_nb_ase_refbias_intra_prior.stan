// negative and beta binomial for ASE eQTL with fixed genotypes but haplotype error accommodating complete allelic imbalance and reference panel bias correction including uncertainty in estimates (version2). Allows for any mixture of priors for eQTL estimate, including single normal
 
data {
  int<lower=0> N; // number  individuals
  int<lower=0> A; // # of individuals with ASE
  int<lower=0> L; // length of vectors with n counts and p(H)
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N]; // total gene counts
  int g[N]; // rnsp geno for all individuals
  int gase[A]; // genotype ASE individuals
  int m[A]; // total ase counts
  int n[L]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  vector[L] ai0; // allelic imbalance estimate for each haplotype in each sample in log scale
  vector[L] sdai0; // standard deviation for allelic imbalance estimate for each haplotype for each sample
  int s[A]; // number of haplotypes per individual
  matrix[N,1+K] cov;
}



parameters {
  vector[K] betas; // regression param
  real bj; // log fold change ASE
  real<lower=1e-5,upper=1e5> phi; //overdipersion param for neg binom
  real<lower=1e-5,upper=1e5> theta; //the overdispersion parameter for beta binomial
  vector[L] rai0; // random intercept AI
}

model {
  // include transformed parameters of no interest
  real ebj;
  vector[k] lps; // help for mixed gaussians
  vector[N] intercept = rep_vector(0, N); // the genetic effect
  real l1pebj;

  // Priors
  theta ~ gamma(1,0.1); //  mean 10 
  phi ~ gamma(1,0.1);  // mean 10

  bj ~ normal(0, 0.05);
  
  // mean expression and covariates
  betas[1] ~ normal(6, 4); // stan normal is mean and sd
  for(i in 2:K) {
    betas[i] ~ cauchy(0, 2.5);//prior for the slopes following Gelman 2008   
  }
  rai0 ~ normal(ai0, sdai0);

  // Likelihood
  ebj = exp(bj); // avoid repeating same calculation
  l1pebj = log1p(ebj) - log(2);

  // Likelihood
  for (i in 1:N) { // log1p(exp(b_j)) - log(2) if het or bj if hom
    if (fabs(g[i]) == 1) {
      intercept[i] = l1pebj;
    }
    if (g[i] == 2) {
      intercept[i] = bj;
    }
  }
  Y ~ neg_binomial_2_log_glm(cov[, 2:cols(cov)], intercept, betas, phi);
}
