functions {
  matrix getSigma(real[] eta_sq, real[] rho_sq, real sigma_sq, matrix[] sqDist, int[] cellLines) {
     matrix[size(cellLines),size(cellLines)] Sigma; 
     int P; 
     int n;
     n <- size(cellLines);
     P <- size(eta_sq);
     for (i in 1:(n-1)) {
        for (j in (i+1):n) {
        real temp[P];
        for (p in 1:P)
          temp[p] <- eta_sq[p] * exp(-rho_sq[p] * sqDist[p][ cellLines[i], cellLines[j] ] );
        Sigma[i,j] <- sum(temp);
        Sigma[j,i] <- Sigma[i,j];
      }
    }
    for (i in 1:n)
      Sigma[i,i] <- sum(eta_sq) + sigma_sq;
    return Sigma; 
  }
}
data {
  int<lower=0> N; // total observations
  int<lower=0> Ntest; 
  int<lower=0> C; 
  int<lower=0> P; 
  vector[N] y;
  int cellLines[N]; 
  int cellLinesTest[Ntest];
  matrix[C,C] sqDist[P]; 
}
parameters {
  real mu; // mean synergy
  real<lower=0> eta_sq[P]; // GP1 signal variance
  real<lower=0> inv_rho_sq[P]; // squared length scale GP 1
  real<lower=0> sigma_sq; // GP 1 noise variance
}
transformed parameters {
  real<lower=0> rho_sq[P]; 
  for (p in 1:P)
    rho_sq[p] <- inv(inv_rho_sq[p]);
}
model {
  matrix[N,N] Sigma; 
  Sigma <- getSigma(eta_sq, rho_sq, sigma_sq, sqDist, cellLines);
  y ~ multi_normal( rep_vector(mu,N) , Sigma);

  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);
}
generated quantities {
  vector[Ntest] ytest; 
  {
    matrix[Ntest,Ntest] Omega;
    matrix[N,Ntest] kappa;
    matrix[Ntest,N] K_transpose_div_Sigma;
    matrix[N,N] Sigma; 

    Sigma <- getSigma( eta_sq, rho_sq, sigma_sq, sqDist, cellLines);
    Omega <- getSigma( eta_sq, rho_sq, sigma_sq, sqDist, cellLinesTest);
    for (i in 1:N)
      for (j in 1:Ntest) {
        real temp[P];
        for (p in 1:P)
          temp[p] <- eta_sq[p] * exp(-rho_sq[p] * sqDist[p][ cellLines[i], cellLinesTest[j] ] );
        kappa[i,j] <- sum(temp);
      }
     K_transpose_div_Sigma <- kappa' / Sigma;
     ytest <- rep_vector(mu,Ntest) + K_transpose_div_Sigma * ( y - rep_vector(mu,N) );
   }
}
