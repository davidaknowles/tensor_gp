data {
  int<lower=0> K; 
  int<lower=0> D; 
  int<lower=0> C;
  int<lower=0> N;
  int<lower=0,upper=C> drugA[N]; 
  int<lower=0,upper=C> drugB[N];
  int<lower=0,upper=C> cellLine[N]; 
  int<lower=0> Ntest;
  int<lower=0,upper=C> drugAtest[Ntest]; 
  int<lower=0,upper=C> drugBtest[Ntest];
  int<lower=0,upper=C> cellLinetest[Ntest]; 
  real synergy[N];
}
parameters {
  real drugV[D,K];
  matrix[C,K] cellV;
  real<lower=0> noise; 
  matrix[D,D] meanSyn;
  real<lower=0> drugscale[K];
  real<lower=0> scalescale;
}
model {
  scalescale ~ cauchy(0,1);
  drugscale ~ cauchy(0,scalescale);
  for (d in 1:D)
    for (k in 1:K)
      drugV[d,k] ~ cauchy(0,drugscale[k]);
  to_vector( cellV ) ~ cauchy(0,1);
  noise ~ cauchy(0,5);
  // to_vector(meanSyn) ~ cauchy(0,5);
  for (n in 1:N) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drugA[n],k] * drugV[drugB[n],k] * cellV[cellLine[n],k]; 
    synergy[n] ~ normal( meanSyn[drugA[n],drugB[n]] + sum(temp), noise );
  }
}
generated quantities {
  real test_synergy[Ntest];
  for (n in 1:Ntest) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drugAtest[n],k] * drugV[drugBtest[n],k] * cellV[cellLinetest[n],k]; 
    test_synergy[n] <- meanSyn[drugA[n],drugB[n]] + sum(temp);
  }
}