getSigma=function(eta_sq_cl, rho_sq_cl, sigma_sq_cl, sqDist_cl, eta_sq_dr, rho_sq_dr, sigma_sq_dr, sqDist_dr, sigma_sq, cellLines, drugA, drugB) {
     n = length(cellLines)
     Sigma=matrix(NA,n,n)
     P=length(eta_sq_cl)
     P <- length(eta_sq_cl);
     P_dr <- length(eta_sq_dr);
     for (i in 1:n) {
        print(i)
        for (j in i:n) {
        temp=numeric(P);
        temp_dr_a=numeric(P_dr)
        temp_dr_b=numeric(P_dr)
        for (p in 1:P)
          temp[p] <- eta_sq_cl[p] * exp(-rho_sq_cl[p] * sqDist_cl[[p]][ cellLines[i], cellLines[j] ] ) ;
        for (p in 1:P_dr)
          temp_dr_a[p] <- eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[[p]][ drugA[i], drugA[j] ] );
        for (p in 1:P_dr)
          temp_dr_b[p] <- eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[[p]][ drugB[i], drugB[j] ] );
        Sigma[i,j] <- ( sum(temp) + ifelse(cellLines[i]==cellLines[j], sigma_sq_cl, 0.0) ) * ( sum(temp_dr_a) + ifelse(drugA[i]==drugA[j], sigma_sq_dr, 0.0) ) * ( sum(temp_dr_b) + ifelse(drugB[i]==drugB[j], sigma_sq_dr, 0.0) ) + ifelse(i==j, sigma_sq, 0.0); 
      }
    }
    for (i in 1:(n-1)) 
        for (j in (i+1):n) 
          Sigma[j,i] <- Sigma[i,j];
    Sigma
}

load("cached_results/sub2_sub2final_tissue1_seed1_iter30.RData")

fit_params=o$par[1:10]

attach(dat)

attach(fit_params)

    kappa=matrix(NA,N,Ntest)
    
    Sigma <- getSigma( eta_sq_cl, rho_sq_cl, sigma_sq_cl, sqDist_cl, eta_sq_dr, rho_sq_dr, sigma_sq_dr, sqDist_dr, sigma_sq, cellLines, drugA, drugB )
    #Omega <- getSigma( eta_sq_cl, rho_sq_cl, sigma_sq_cl, sqDist_cl, eta_sq_dr, rho_sq_dr, sigma_sq_dr, sqDist_dr, sigma_sq, cellLinesTest, drugATest, drugBTest)
    for (i in 1:N)
      for (j in 1:Ntest) {
        temp=numeric(P)
        temp_dr_a=numeric(P_dr)
  		  temp_dr_b=numeric(P_dr)
        for (p in 1:P)
          temp[p] <- eta_sq_cl[p] * exp(-rho_sq_cl[p] * sqDist_cl[[p]][ cellLines[i], cellLinesTest[j] ] ) ;
        for (p in 1:P_dr)
          temp_dr_a[p] <- eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[[p]][ drugA[i], drugATest[j] ] );
        for (p in 1:P_dr)
          temp_dr_b[p] <- eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[[p]][ drugB[i], drugBTest[j] ] );
        kappa[i,j] <- ( sum(temp) + ifelse(cellLines[i]==cellLinesTest[j], sigma_sq_cl, 0.0) ) * ( sum(temp_dr_a) + ifelse(drugA[i]==drugATest[j], sigma_sq_dr, 0.0) ) * ( sum(temp_dr_b) + ifelse(drugB[i]==drugBTest[j], sigma_sq_dr, 0.0) ) ; 
      }
     K_transpose_div_Sigma <- solve(Sigma, kappa)
     ytest <- rep(mu,Ntest) + K_transpose_div_Sigma * ( y - rep(mu,N) );

