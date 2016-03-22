
getSigmaFast=function(eta_sq_cl, rho_sq_cl, sigma_sq_cl, sqDist_cl, eta_sq_dr, rho_sq_dr, sigma_sq_dr, sqDist_dr, sigma_sq, cellLines, drugA, drugB) {
  P <- length(eta_sq_cl);
  P_dr <- length(eta_sq_dr);
  temp=Reduce("+", foreach(p=1:P) %dopar% { eta_sq_cl[p] * exp(-rho_sq_cl[p] * sqDist_cl[[p]][ cellLines, cellLines ] ) } )
  temp_dr_a =Reduce("+", foreach(p=1:P_dr) %dopar% { eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[[p]][ drugA, drugA ] ) } )
  temp_dr_b =Reduce("+", foreach(p=1:P_dr) %dopar% { eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[[p]][ drugB, drugB ] ) } )
  (temp + sigma_sq_cl*outer(cellLines,cellLines,"==")) * (temp_dr_a + sigma_sq_dr*outer(drugA,drugA,"==")) * (temp_dr_b + sigma_sq_dr*outer(drugB,drugB,"==")) + sigma_sq * diag(length(cellLines)) 
}

getKappaFast=function(eta_sq_cl, rho_sq_cl, sigma_sq_cl, sqDist_cl, eta_sq_dr, rho_sq_dr, sigma_sq_dr, sqDist_dr, sigma_sq, cellLines, drugA, drugB, cellLinesTest, drugATest, drugBTest) {
  P <- length(eta_sq_cl);
  P_dr <- length(eta_sq_dr);
  temp=Reduce("+", foreach(p=1:P) %dopar% { eta_sq_cl[p] * exp(-rho_sq_cl[p] * sqDist_cl[[p]][ cellLines, cellLinesTest ] ) } )
  temp_dr_a =Reduce("+", foreach(p=1:P_dr) %dopar% { eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[[p]][ drugA, drugATest ] ) } )
  temp_dr_b =Reduce("+", foreach(p=1:P_dr) %dopar% { eta_sq_dr[p] * exp(-rho_sq_dr[p] * sqDist_dr[[p]][ drugB, drugBTest ] ) } )
  (temp + sigma_sq_cl*outer( cellLines, cellLinesTest, "==" ) ) * (temp_dr_a + sigma_sq_dr*outer( drugA, drugATest, "==" ) ) * (temp_dr_b + sigma_sq_dr*outer( drugB, drugBTest, "==" )) 
}
