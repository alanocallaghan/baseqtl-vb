// negative and beta binomial for ASE eQTL with fixed genotypes but haplotype error accommodating complete allelic imbalance, version 2. Allows for any mixture of gaussians for bj prior (eQTL effect).
 
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
  int s[A]; // number of haplotypes per individual
  matrix[N,1+K] cov;
  vector[k] aveP; // mean for prior Gaussians for eQTL effect prior
  vector[k] sdP; // sd for prior Gaussians for eQTL effect prior
  vector[k] mixP; // log of mixing proportions for eQTL effect prior

  }

parameters {
  vector[K] betas; // regression param
  real bj; // log fold change ASE
  real<lower=0> phi; //overdipersion param for neg binom
  real<lower=0> theta; //the overdispersion parameter for beta binomial
}

model {
  // include transformed parameters of no interest
  vector[N] lmu;//the linear predictor
  vector[A] p; // ASE proportion
  real ebj;
  real debj;
  int pos; // to loop over haplotypes for each individual
  vector[L] ltmp; //  log BB likelihood
  vector[k] lps; // help for mixed gaussians


  // Priors
  theta ~ gamma(1,0.1); //  mean 10 
  phi ~ gamma(1,0.1);  // mean 10
  betas[1] ~ normal(6,4); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008   
      }
  // mixture of gaussians for bj:
  for(i in 1:k){
    lps[i] = normal_lpdf(bj | aveP[i], sdP[i]) + mixP[i];
  }
  target += log_sum_exp(lps);


  // Likelihood
  ebj=exp(bj); // avoid repeating same calculation
  debj=ebj/(1+ebj);

  lmu = cov[,2:cols(cov)]*betas;
  for(i in 1:N){ // neg binomial
    lmu[i] = fabs(g[i])==1 ? lmu[i] + log1p(ebj)-log(2) : lmu[i];
    lmu[i] = g[i]==2 ? lmu[i] + bj : lmu[i];
    target += neg_binomial_2_lpmf(Y[i] | exp(lmu[i]),phi);
  }
  
  pos = 1;
  for(i in 1:A){ // ASE
    p[i]= gase[i]==1 ? debj : 0.5;
    p[i]= gase[i]==-1 ? 1-debj : p[i];  // haplotype swap
     for (r in pos:(pos+s[i]-1)){
       ltmp[r]=beta_binomial_lpmf(n[r] | m[i], p[i]*theta , (1-p[i])*theta) + log(pH[r]);
     }
     target += log_sum_exp(ltmp[pos:(pos+s[i]-1)]);
     pos=pos+s[i];
 
  }
}
