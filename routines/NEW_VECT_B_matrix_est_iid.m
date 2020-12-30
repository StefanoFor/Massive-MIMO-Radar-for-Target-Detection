function [B_est] = NEW_VECT_B_matrix_est_iid(N,v,x,alpha,g_f_R_mat,g_f_I_mat,ii,jj,I,J)

d = 2;

hat_c = x - alpha*v;

c_R = real(hat_c);
c_I = imag(hat_c);

C_RR = [c_R';c_R'].^2;
C_II = [c_R';c_R'].^2;

Sigma_R_full = c_R*c_R';
Sigma_I_full = c_I*c_I';

Sigma_R = Sigma_R_full - Sigma_R_full.*(ii - I(1) > (diff(I(1:2)))/(diff(J(1:2)))*jj | ii - I(end) + 1< (diff(I(3:end)))/(diff(J(3:end)))*(jj - J(end)+1));
B_2_R = 4*N^(-1)*g_f_R_mat*Sigma_R*g_f_R_mat';

Sigma_I = Sigma_I_full - Sigma_I_full.*(ii - I(1) > (diff(I(1:2)))/(diff(J(1:2)))*jj | ii - I(end) + 1< (diff(I(3:end)))/(diff(J(3:end)))*(jj - J(end)+1));
B_2_I = 4*N^(-1)*g_f_I_mat*Sigma_I*g_f_I_mat';

clear Sigma_R_full Sigma_I_full Sigma_R Sigma_I

B_1_R = 4*N^(-1)*(C_RR.*g_f_R_mat)*g_f_R_mat';
B_1_I = 4*N^(-1)*(C_II.*g_f_I_mat)*g_f_I_mat';

% B_2_R = zeros(d,d);
% B_2_I = zeros(d,d);
% for m = 1:M
%     for n = m+1 : N
%         B_2_R = B_2_R + 4*N^(-1)* c_R(n) * c_R(n-m) * ( g_f_R_mat(:,n)*g_f_R_mat(:,n-m)' + g_f_R_mat(:,n-m)*g_f_R_mat(:,n)' );
%         B_2_I = B_2_I + 4*N^(-1)* c_I(n) * c_I(n-m) * ( g_f_I_mat(:,n)*g_f_I_mat(:,n-m)' + g_f_I_mat(:,n-m)*g_f_I_mat(:,n)' );
%     end
% end

% for m = 1:M
%     for n = m+1 : 2*N
%         if (n-m) <= N && n<= N
%             B_2 = B_2 + 4*N^(-1)* aug_c(n) * aug_c(n-m) * ( g_f_n(n,v)*g_f_n(n-m,v)' + g_f_n(n-m,v)*g_f_n(n,v)' );
%         end
%         if (n-m) > N && n<= N
%             B_2 = B_2 + 4*N^(-1)* aug_c(n) * aug_c(n-m) * ( g_f_n(n,v)*g_f_n(n-m,v)' + g_f_n(n-m,v)*g_f_n(n,v)' );
%         end
%         if (n-m) <= N && n > N
%             B_2 = B_2 + 4*N^(-1)* aug_c(n) * aug_c(n-m) * ( g_f_n(n,v)*g_f_n(n-m,v)' + g_f_n(n-m,v)*g_f_n(n,v)' );
%         end
%         if (n-m) > N && n > N
%             B_2 = B_2 + 4*N^(-1)* aug_c(n) * aug_c(n-m) * ( g_f_n(n,v)*g_f_n(n-m,v)' + g_f_n(n-m,v)*g_f_n(n,v)' );
%         end
%     end
% end


B_est = B_1_R + B_1_I + B_2_R + B_2_I;
end
