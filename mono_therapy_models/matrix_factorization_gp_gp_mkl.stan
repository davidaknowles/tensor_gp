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
  int<lower=0> P; // num views
  int<lower=0> K; // num latent characteristics
  int<lower=0> D; // num drugs
  int<lower=0> C; // num cell lines
  int<lower=0> N; // num measurements
  int<lower=0,upper=C> drug[N]; // indices of drugs
  int<lower=0,upper=C> cellLine[N]; // indices of cell lines
  vector[N] y; // sensitivity measurements (e.g. Einf or IC50)
  int<lower=0> Ntest; // number of test measurements
  int<lower=0,upper=C> drugtest[Ntest]; // test indices
  int<lower=0,upper=C> cellLinetest[Ntest];
  matrix[C,C] sqDist[P]; // squared distance between cell lines (e.g. in gene expression space)
}
parameters {
  real drugV[D,K]; // latent values for drugs
  vector[C] z[K]; // whitened latent values for cell lines
  real meanSyn; // mean sensitivity
  real<lower=0> drugscale[K]; // scale for each factor (ARD)
  real<lower=0> scalescale; // overall scale
  real<lower=0> eta_sq[P]; // GP1 signal variance
  real<lower=0> inv_rho_sq[P]; // squared length scale GP 1
  real<lower=0> sigma_sq; // GP 1 noise variance
  real<lower=0> eta_sq_tau; // GP 2 signal variance
  real<lower=0> inv_rho_sq_tau; // GP 2 squared length scale
  real<lower=0> sigma_sq_tau; // GP 2 noise variance
}
transformed parameters {
  vector[C] cellV[K]; // cell line latent values (derived from z)
  matrix[C,C] Sigma; // GP1 Gram matrix
  matrix[C,C] L; // Cholesky of Sigma
  real<lower=0> rho_sq[P]; // GP1 inverse squared length scale
  real<lower=0> rho_sq_tau; // GP2 inverse squared length scale
  vector[2*K] x[N]; // concatenation of cellV and drugV for each measurement
  vector[2*K] xtest[Ntest]; // concatenation of cellV and drugV for each test measurement
  vector[N] mu; // repeated "mean"
  
  rho_sq_tau <- inv(inv_rho_sq_tau);

  // build GP1's Gram matrix
  for (p in 1:P)
    rho_sq[p] <- inv(inv_rho_sq[p]);
  for (i in 1:(C-1)) {
    for (j in (i+1):C) {
      real temp[P];
      for (p in 1:P)
        temp[p] <- eta_sq[p] * exp(-rho_sq[p] * sqDist[p][i,j]);
      Sigma[i,j] <- sum(temp);
      Sigma[j,i] <- Sigma[i,j];
    } 
  }
  for (k in 1:C)
    Sigma[k,k] <- sum(eta_sq) + sigma_sq;
  L <- cholesky_decompose(Sigma);
  for (k in 1:K) 
    cellV[k] <- L * z[k];
  for (n in 1:N) {
    for (k in 1:K) {
      x[n][k] <- drugV[drug[n],k];
      x[n][K+k] <- cellV[k][cellLine[n] ];
    }
  }
  for (n in 1:N) 
    mu[n] <- meanSyn; 
  for (n in 1:Ntest) {
    for (k in 1:K) {
      xtest[n][k] <- drugV[drugtest[n],k];
      xtest[n][K+k] <- cellV[k][cellLinetest[n] ];
    }
  }
}
model {
  // Build GP2's Gram matrix
  matrix[N,N] Tau;
  Tau  <- getTau(eta_sq_tau, rho_sq_tau, sigma_sq_tau, x);

  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);

  eta_sq_tau ~ cauchy(0,5);
  inv_rho_sq_tau ~ cauchy(0,5);
  sigma_sq_tau ~ cauchy(0,5);

  scalescale ~ gamma(1,1);
  drugscale ~ gamma(scalescale/K,scalescale);
  for (d in 1:D)
    for (k in 1:K)
      drugV[d,k] ~ normal(0,drugscale[k]);
  
  y ~ multi_normal( mu , Tau);
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
    Tau <- getTau(eta_sq_tau, rho_sq_tau, sigma_sq_tau, x);
    Omega <- getTau( eta_sq_tau, rho_sq_tau, sigma_sq_tau, xtest);

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