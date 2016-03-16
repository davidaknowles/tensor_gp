functions {
  matrix getTau(real eta_sq_tau, real rho_sq_tau, real sigma_sq_tau, vector[] x) {
    int N;
    matrix[size(x),size(x)] Tau;
    N <- size(x);
    
    for (n in 1:(N-1)) {
      for (m in (n+1):N) {
        Tau[m,n] <- eta_sq_tau * exp(-rho_sq_tau * dot_self(x[m] - x[n]));
        Tau[n,m] <- Tau[m,n];
      }
    }
    for (n in 1:N) 
      Tau[n,n] <- eta_sq_tau + sigma_sq_tau;
    return Tau; 
  }
}
data {
  int<lower=0> K; // num latent characteristics
  int<lower=0> D; // num drugs
  int<lower=0> C; // num cell lines
  int<lower=0> Nmono; // num measurements
  int<lower=0,upper=C> drugMono[Nmono]; // indices of drugs
  int<lower=0,upper=C> cellLineMono[Nmono]; // indices of cell lines
  vector[Nmono] y; // monotherapy sensitivity measurements (e.g. Einf or IC50)
  
  int<lower=0> N; // num measurements
  int<lower=0,upper=C> drugA[N]; // indices of drug A
  int<lower=0,upper=C> drugB[N]; // indices of drug B
  int<lower=0,upper=C> cellLine[N]; // indices of cell lines
  vector[N] syn; // synergy measurements

  int<lower=0> Ntest; // number of test measurements
  int<lower=0,upper=C> drugAtest[Ntest]; // test indices
  int<lower=0,upper=C> drugBtest[Ntest];
  int<lower=0,upper=C> cellLinetest[Ntest];

  matrix[C,C] sqDist; // squared distance between cell lines (e.g. in gene expression space)
}
parameters {
  real drugV[D,K]; // latent values for drugs
  vector[C] z[K]; // whitened latent values for cell lines

  real meanSens; // mean sensitivity
  real meanSyn; // mean synergy
  real<lower=0> drugscale[K]; // scale for each factor (ARD)
  real<lower=0> scalescale; // overall scale

  real<lower=0> eta_sq; // GP1 signal variance
  real<lower=0> inv_rho_sq; // squared length scale GP 1
  real<lower=0> sigma_sq; // GP 1 noise variance

  real<lower=0> eta_sq_tau; // GP 2 signal variance
  real<lower=0> inv_rho_sq_tau; // GP 2 squared length scale
  real<lower=0> sigma_sq_tau; // GP 2 noise variance

  real<lower=0> eta_sq_delta; // GP 3 signal variance
  real<lower=0> inv_rho_sq_delta; // GP 3 squared length scale
  real<lower=0> sigma_sq_delta; // GP 3 noise variance
}
transformed parameters {
  vector[C] cellV[K]; // cell line latent values (derived from z)
  matrix[C,C] Sigma; // GP1 Gram matrix
  matrix[C,C] L; // Cholesky of Sigma
  real<lower=0> rho_sq; // GP1 inverse squared length scale
  real<lower=0> rho_sq_tau; // GP2 inverse squared length scale
  real<lower=0> rho_sq_delta; // GP3 inverse squared length scale
  vector[3*K] x[N]; // concatenation of cellV and drugV for each measurement
  vector[3*K] xtest[Ntest]; // concatenation of cellV and drugV for each test measurement
  vector[2*K] xmono[Nmono];
  vector[Nmono] muMono; // repeated "mean"
  vector[N] muSyn; 

  rho_sq_tau <- inv(inv_rho_sq_tau);
  rho_sq_delta <- inv(inv_rho_sq_delta);

  // build GP1's Gram matrix
  rho_sq <- inv(inv_rho_sq);
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
      x[n][k] <- drugV[drugA[n],k];
      x[n][K+k] <- drugV[drugB[n],k];
      x[n][2*K+k] <- cellV[k][cellLine[n] ];
    }
  }
  for (n in 1:N) 
    muSyn[n] <- meanSyn; 
  for (n in 1:Nmono) 
    muMono[n] <- meanSens; 
  for (n in 1:Ntest) {
    for (k in 1:K) {
      xtest[n][k] <- drugV[drugAtest[n],k];
      xtest[n][K+k] <- drugV[drugBtest[n],k];
      xtest[n][2*K+k] <- cellV[k][cellLinetest[n] ];
    }
  }
  for (n in 1:Nmono) {
    for (k in 1:K) {
      xmono[n][k] <- drugV[drugMono[n],k];
      xmono[n][K+k] <- cellV[k][cellLineMono[n] ];
    }
  }
}
model {
  // Build GP2's Gram matrix
  matrix[Nmono,Nmono] Tau;
  matrix[N, N] Delta;   

  Tau  <- getTau(eta_sq_tau, rho_sq_tau, sigma_sq_tau, xmono);

  Delta <- getTau(eta_sq_delta, rho_sq_delta, sigma_sq_delta, x);
  
  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);

  eta_sq_tau ~ cauchy(0,5);
  inv_rho_sq_tau ~ cauchy(0,5);
  sigma_sq_tau ~ cauchy(0,5);

  eta_sq_delta ~ cauchy(0,5);
  inv_rho_sq_delta ~ cauchy(0,5);
  sigma_sq_delta ~ cauchy(0,5);

  scalescale ~ cauchy(0,1);
  drugscale ~ cauchy(0,scalescale);
  for (d in 1:D)
    for (k in 1:K)
      drugV[d,k] ~ cauchy(0,drugscale[k]);

  y ~ multi_normal( muMono, Tau);
  syn ~ multi_normal( muSyn, Delta); 
}
generated quantities {
  // get predictions on test measurements
  vector[Ntest] testMean;
  vector[Ntest] testVar;
  {
    matrix[Ntest,Ntest] testVarMat;
    matrix[Ntest,Ntest] Omega;
    matrix[N,Ntest] Kmat;
    matrix[Ntest,N] K_transpose_div_Sigma;
    matrix[N,N] Tau;
    Tau <- getTau(eta_sq_delta, rho_sq_delta, sigma_sq_delta, x);
    Omega <- getTau( eta_sq_delta, rho_sq_delta, sigma_sq_delta, xtest);

    for (n in 1:N) 
      for (m in 1:Ntest) 
        Kmat[n,m] <- eta_sq_tau * exp(-rho_sq_tau * dot_self(x[n] - xtest[m]));
    K_transpose_div_Sigma <- Kmat' / Tau;
    // vector[Ntest] muTest; 
    testMean <- meanSyn + K_transpose_div_Sigma * (syn - muSyn);
    testVarMat <- Omega - K_transpose_div_Sigma * Kmat;
    // for (i in 1:(Ntest-1))
    //  for (j in (i+1):Ntest)
    //    testVar[i,j] <- testVar[j,i];
    testVar <- diagonal(testVarMat);
  }
}