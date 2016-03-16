functions {
  matrix getSigma(real[] eta_sq_cl, real[] rho_sq_cl, real sigma_sq_cl, matrix[] sqDist_cl, real[] eta_sq_dr, real[] rho_sq_dr, real sigma_sq_dr, matrix[] sqDist_dr, real sigma_sq, int[] cellLines, int[] drugs) {
     matrix[size(cellLines),size(cellLines)] Sigma; 
     int P; 
     int P_dr;
     int n;
     n <- size(cellLines);
     P <- size(eta_sq_cl);
     P_dr <- size(eta_sq_dr);
     for (i in 1:n) {
        for (j in i:n) {
        real temp[P];
        real temp_dr[P_dr];
        for (p in 1:P)
          temp[p] <- eta_sq_cl[p] * exp(-rho_sq_cl[p] * sqDist_cl[p][ cellLines[i], cellLines[j] ] ) ;
        for (p in 1:P_dr)
          temp_dr[p] <- eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[p][ drugs[i], drugs[j] ] );
        Sigma[i,j] <- ( sum(temp) + if_else(cellLines[i]==cellLines[j], sigma_sq_cl, 0.0) ) * ( sum(temp_dr) + if_else(drugs[i]==drugs[j], sigma_sq_dr, 0.0) ) + if_else(i==j, sigma_sq, 0.0); 
        Sigma[j,i] <- Sigma[i,j];
      }
    }
    return Sigma; 
  }
}
data {
  int<lower=0> N; // total observations
  int<lower=0> Ntest; 
  int<lower=0> C; 
  int<lower=0> D; 
  int<lower=0> P; 
  int<lower=0> P_dr; 
  vector[N] y;
  int cellLines[N]; 
  int cellLinesTest[Ntest];
  int drugs[N]; 
  int drugsTest[Ntest];
  matrix[C,C] sqDist_cl[P]; 
  matrix[D,D] sqDist_dr[P_dr]; 
}
parameters {
  real mu; // mean synergy
  real<lower=0> eta_sq_cl[P]; // GP1 signal variance
  real<lower=0> inv_rho_sq_cl[P]; // squared length scale GP 1
  real<lower=0> sigma_sq_cl; // GP 1 noise variance
  real<lower=0> eta_sq_dr[P_dr]; // GP1 signal variance
  real<lower=0> inv_rho_sq_dr[P_dr]; // squared length scale GP 1
  real<lower=0> sigma_sq_dr; // GP 1 noise variance
  real<lower=0> sigma_sq; 
}
transformed parameters {
  real<lower=0> rho_sq_dr[P_dr]; 
  real<lower=0> rho_sq_cl[P]; 
  for (p in 1:P)
    rho_sq_cl[p] <- inv(inv_rho_sq_cl[p]);
  for (p in 1:P_dr)
    rho_sq_dr[p] <- inv(inv_rho_sq_dr[p]);
}
model {
  matrix[N,N] Sigma; 
  Sigma <- getSigma(eta_sq_cl, rho_sq_cl, sigma_sq_cl, sqDist_cl, eta_sq_dr, rho_sq_dr, sigma_sq_dr, sqDist_dr, sigma_sq, cellLines, drugs);
  y ~ multi_normal( rep_vector(mu,N) , Sigma);

  eta_sq_cl ~ cauchy(0,5);
  inv_rho_sq_cl ~ cauchy(0,5);
  sigma_sq_cl ~ cauchy(0,5);
  
  eta_sq_dr ~ cauchy(0,5);
  inv_rho_sq_dr ~ cauchy(0,5);
  sigma_sq_dr ~ cauchy(0,5);
  
  sigma_sq ~ cauchy(0,5);
}
generated quantities {
  vector[Ntest] ytest; 
  {
    matrix[Ntest,Ntest] Omega;
    matrix[N,Ntest] kappa;
    matrix[Ntest,N] K_transpose_div_Sigma;
    matrix[N,N] Sigma; 

    Sigma <- getSigma( eta_sq_cl, rho_sq_cl, sigma_sq_cl, sqDist_cl, eta_sq_dr, rho_sq_dr, sigma_sq_dr, sqDist_dr, sigma_sq, cellLines, drugs );
    Omega <- getSigma( eta_sq_cl, rho_sq_cl, sigma_sq_cl, sqDist_cl, eta_sq_dr, rho_sq_dr, sigma_sq_dr, sqDist_dr, sigma_sq, cellLinesTest, drugsTest);
    for (i in 1:N)
      for (j in 1:Ntest) {
        real temp[P];
        real temp_dr[P_dr];
        for (p in 1:P)
          temp[p] <- eta_sq_cl[p] * exp(-rho_sq_cl[p] * sqDist_cl[p][ cellLines[i], cellLinesTest[j] ] ) ;
        for (p in 1:P_dr)
          temp_dr[p] <- eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[p][ drugs[i], drugsTest[j] ] );
        kappa[i,j] <- ( sum(temp) + if_else(cellLines[i]==cellLinesTest[j], sigma_sq_cl, 0.0) ) * ( sum(temp_dr) + if_else(drugs[i]==drugsTest[j], sigma_sq_dr, 0.0) ) ; 
      }
     K_transpose_div_Sigma <- kappa' / Sigma;
     ytest <- rep_vector(mu,Ntest) + K_transpose_div_Sigma * ( y - rep_vector(mu,N) );
   }
}
