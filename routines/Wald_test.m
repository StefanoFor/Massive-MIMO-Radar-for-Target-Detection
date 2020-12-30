function W_T = Wald_test(N,v,x,hat_alpha,l,n2)

B_est_mod = B_matrix_est(N,v,x,hat_alpha,l);
    
% Wald Test

W_T = 2 * n2^2 * (abs(hat_alpha)^2) /B_est_mod;

end

