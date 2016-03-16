data {
  int<lower=0> K; 
  int<lower=0> D; 
  int<lower=0> C;
  int<lower=0> N;
  int<lower=0,upper=C> drug[N]; 
  int<lower=0,upper=C> cellLine[N]; 
  int<lower=0> Ntest;
  int<lower=0,upper=C> drugtest[Ntest]; 
  int<lower=0,upper=C> cellLinetest[Ntest];
  matrix[C,C] sqDist; 
  vector[N] y;
}
parameters {
  real drugV[D,K];
  vector[C] z[K];
  real<lower=0> noise; 
  real meanSyn;
  real<lower=0> drugscale[K];
  real<lower=0> scalescale;
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
  real<lower=0> eta_sq_tau;
  real<lower=0> inv_rho_sq_tau;
  
  real<lower=0> sigma_sq_tau;
}
transformed parameters {
  matrix[N,N] Tau;
  vector[C] cellV[K];
  matrix[C,C] Sigma;
  matrix[C,C] L;
  real<lower=0> rho_sq;
  real<lower=0> rho_sq_tau;
  vector[2*K] x[N];
  vector[2*K] xtest[Ntest];
  vector[N] mu; 
  rho_sq <- inv(inv_rho_sq);
  rho_sq_tau <- inv(inv_rho_sq_tau);
  for (i in 1:(C-1)) {
    for (j in (i+1):C) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * sqDist[i,j]);
      Sigma[j,i] <- Sigma[i,j];
    } 
  }
  for (k in 1:C)
    Sigma[k,k] <- eta_sq + sigma_sq;
  L <- cholesky_decompose(Sigma);
  for (k in 1:K) 
    cellV[k] <- L * z[k];
  for (n in 1:N) {
    for (k in 1:K) {
      x[n][k] <- drugV[drug[n],k];
      x[n][K+k] <- cellV[k][cellLine[n] ];
    }
  }
  for (n in 1:Ntest) {
    for (k in 1:K) {
      xtest[n][k] <- drugV[drugtest[n],k];
      xtest[n][K+k] <- cellV[k][cellLinetest[n] ];
    }
  }
  for (n in 1:(N-1)) {
    for (m in (n+1):N) {
      Tau[m,n] <- eta_sq_tau * exp(-rho_sq_tau * dot_self(x[m] - x[n]));
      Tau[n,m] <- Tau[m,n];
    }
  }
  for (n in 1:N) {
    Tau[n,n] <- eta_sq_tau + sigma_sq_tau;
    mu[n] <- meanSyn;
  }
}
model {
  
  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);

  eta_sq_tau ~ cauchy(0,5);
  inv_rho_sq_tau ~ cauchy(0,5);
  sigma_sq_tau ~ cauchy(0,5);

  scalescale ~ cauchy(0,1);
  drugscale ~ cauchy(0,scalescale);
  for (d in 1:D)
    for (k in 1:K)
      drugV[d,k] ~ cauchy(0,drugscale[k]);
  
  noise ~ cauchy(0,5);

  y ~ multi_normal( mu , Tau);
}
generated quantities {
  vector[Ntest] testMean;
  vector[Ntest] testVar;
  {
    matrix[Ntest,Ntest] testVarMat;
    matrix[Ntest,Ntest] Omega;
    matrix[N,Ntest] Kmat;
    matrix[Ntest,N] K_transpose_div_Sigma;

    for (n in 1:(Ntest-1)) 
      for (m in (n+1):Ntest) {
        Omega[m,n] <- eta_sq_tau * exp(-rho_sq_tau * dot_self(xtest[m] - xtest[n]));
        Omega[n,m] <- Omega[m,n];
      }
    for (n in 1:Ntest) 
      Omega[n,n] <- eta_sq_tau + sigma_sq_tau;
    for (n in 1:N) 
      for (m in 1:Ntest) 
        Kmat[n,m] <- eta_sq_tau * exp(-rho_sq_tau * dot_self(x[n] - xtest[m]));
    K_transpose_div_Sigma <- Kmat' / Tau;
    // vector[Ntest] muTest; 
    testMean <- meanSyn + K_transpose_div_Sigma * (y - mu);
    testVarMat <- Omega - K_transpose_div_Sigma * Kmat;
    // for (i in 1:(Ntest-1))
    //  for (j in (i+1):Ntest)
    //    testVar[i,j] <- testVar[j,i];
    testVar <- diagonal(testVarMat);
  }
}