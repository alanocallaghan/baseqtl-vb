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
  vector[L] p; // ASE proportion
  int pos; // to loop over haplotypes for each individual
  vector[L] ltmp; //  log BB likelihood
  vector[L] esum; // reduce computation inverse logit (rai0 + bj)
  vector[L] esum0; // allelic imbalance proportion under the null
  vector[k] lps; // help for mixed gaussians

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
  pos = 1;
  esum = inv_logit(rai0 + bj);
  esum0 = inv_logit(rai0);
  
  for(i in 1:A) { // ASE
    for (r in pos:(pos+s[i]-1)) {
       
      p[r]= gase[i]==1 ? esum[r] : esum0[r];
      p[r]= gase[i]==-1 ? 1-esum[r] : p[r];  // haplotype swap
       
      ltmp[r]=beta_binomial_lpmf(n[r] | m[i], p[r]*theta , (1-p[r])*theta) + log(pH[r]);
    }
    target += log_sum_exp(ltmp[pos:(pos+s[i]-1)]);
    pos=pos + s[i];
  }
}
