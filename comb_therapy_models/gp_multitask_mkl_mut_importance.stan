
data {
  int<lower=0> N; // total observations
  int<lower=0> num_genes; 
  int<lower=0> C; // number of cell lines
  int<lower=0> D; // number of drugs
  int<lower=0> P; 
  int<lower=0> P_dr; 
  int cellLines[N]; 
  int cellLinesTest;
  int drugA[N]; 
  int drugATest;
  int drugB[N]; 
  int drugBTest;
  matrix[C,C] sqDist_cl[P]; 
  matrix[D,D] sqDist_dr[P_dr]; 
  real<lower=0> eta_sq_cl[P]; // GP1 signal variance
  real<lower=0> inv_rho_sq_cl[P]; // squared length scale GP 1
  real<lower=0> sigma_sq_cl; // GP 1 noise variance
  real<lower=0> eta_sq_cl_mut; // GP1 signal variance
  real<lower=0> inv_rho_sq_mut; // squared length scale GP 1
  real<lower=0> eta_sq_dr[P_dr]; // GP1 signal variance
  real<lower=0> inv_rho_sq_dr[P_dr]; // squared length scale GP 1
  real<lower=0> sigma_sq_dr; // GP 1 noise variance
  vector[num_genes] mut[C];
  vector[N] unmixed_y; 
}
parameters {
  vector[num_genes] mut_test;
}
transformed parameters {
  real<lower=0> rho_sq_dr[P_dr]; 
  real<lower=0> rho_sq_cl[P]; 
  real<lower=0> rho_sq_mut; 
  for (p in 1:P)
    rho_sq_cl[p] <- inv(inv_rho_sq_cl[p]);
  for (p in 1:P_dr)
    rho_sq_dr[p] <- inv(inv_rho_sq_dr[p]);
  rho_sq_mut <- inv_rho_sq_mut;
}
model {
  vector[N] kappa;
  for (i in 1:N) {
      real temp[P+1];
      real temp_dr_a[P_dr];
	    real temp_dr_b[P_dr];
	    real eucl_dist2; 
      for (p in 1:P)
        temp[p] <- eta_sq_cl[p] * exp(-rho_sq_cl[p] * sqDist_cl[p][ cellLines[i], cellLinesTest ] ) ;
      eucl_dist2 <- dot_self( mut[ cellLines[i] ] - mut_test ) ; 
      temp[P+1] <- eta_sq_cl_mut * exp(-rho_sq_mut * (eucl_dist2 / (eucl_dist2 + dot_product(mut[ cellLines[i] ], mut_test )))^2);
      for (p in 1:P_dr)
        temp_dr_a[p] <- eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[p][ drugA[i], drugATest ] );
      for (p in 1:P_dr)
        temp_dr_b[p] <- eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[p][ drugB[i], drugBTest ] );
      kappa[i] <- ( sum(temp) + if_else(cellLines[i]==cellLinesTest, sigma_sq_cl, 0.0) ) * ( sum(temp_dr_a) + if_else(drugA[i]==drugATest, sigma_sq_dr, 0.0) ) * ( sum(temp_dr_b) + if_else(drugB[i]==drugBTest, sigma_sq_dr, 0.0) ) ; 
    }
  increment_log_prob( dot_product(unmixed_y, kappa) );
}
