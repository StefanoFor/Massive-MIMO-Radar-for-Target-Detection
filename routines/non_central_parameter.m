function W_T = non_central_parameter(N,v,x_H0,alpha,g_f_R_mat,g_f_I_mat,l,n2)

B_est_mod = B_matrix_est_iid(N,v,x_H0,0,g_f_R_mat,g_f_I_mat,l);
    
% Wald Test 
inv_B_mod = (B_est_mod(1,1)*B_est_mod(2,2)-B_est_mod(1,2)*B_est_mod(2,1))^(-1)*[B_est_mod(2,2) -B_est_mod(1,2);
                                                                                -B_est_mod(2,1) B_est_mod(1,1)];    
 
theta = [real(alpha); imag(alpha)];

W_T = n2^2 * theta.' * inv_B_mod * theta;