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
  int<lower=0> K; // groups (# drug combinations)
  int<lower=0> C; 
  int<lower=0> P; 
  vector[N] y;
  int cellLines[N]; 
  int cellLinesTest[Ntest];
  int s[K];
  int stest[K];
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
  int pos; 
  pos <- 1; 
  for (k in 1:K) {
     matrix[s[k],s[k]] Sigma; 
     int cellLinesHere[s[k]]; 
     for (i in 1:s[k])
       cellLinesHere[i] <- cellLines[pos+i-1];
     Sigma <- getSigma(eta_sq, rho_sq, sigma_sq, sqDist, cellLinesHere);
     segment(y, pos, s[k]) ~ multi_normal( rep_vector(mu,s[k]) , Sigma);
     pos <- pos + s[k];
   }
  
  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);
}
generated quantities {
  int pos; 
  int postest; 
  vector[Ntest] ytest; 
  pos <- 1; 
  postest <- 1; 
  for (k in 1:K) {
     int n; 
     int ntest;
     vector[stest[k]] pred; 
      matrix[stest[k],stest[k]] Omega;
      matrix[s[k],stest[k]] kappa;
      matrix[stest[k],s[k]] K_transpose_div_Sigma;
      matrix[s[k],s[k]] Sigma; 
     int cellLinesHere[s[k]]; 
     int cellLinesTestHere[stest[k]];

     n <- s[k];
     ntest <- stest[k];

     for (i in 1:n)
       cellLinesHere[i] <- cellLines[pos+i-1];
     for (i in 1:ntest)
       cellLinesTestHere[i] <- cellLines[postest+i-1];
     
     Sigma <- getSigma( eta_sq, rho_sq, sigma_sq, sqDist, cellLinesHere);
     Omega <- getSigma( eta_sq, rho_sq, sigma_sq, sqDist, cellLinesTestHere);
    for (i in 1:n)
      for (j in 1:ntest) {
        real temp[P];
        for (p in 1:P)
          temp[p] <- eta_sq[p] * exp(-rho_sq[p] * sqDist[p][ cellLines[pos+i-1], cellLinesTest[postest+j-1] ] );
        kappa[i,j] <- sum(temp);
      }
     K_transpose_div_Sigma <- kappa' / Sigma;
     pred <- rep_vector(mu,ntest) + K_transpose_div_Sigma * ( segment(y, pos, n) - rep_vector(mu,s[k]) );
     for (j in 1:ntest)
       ytest[postest+j-1] <- pred[j];
     pos <- pos + n;
     postest <- postest + ntest;
   }
}
