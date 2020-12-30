function [B_est] = B_matrix_est(N,v,x,alpha,l)

hat_c = x - alpha*v;

B_1 = sum( (abs(v).^2) .* (abs(hat_c).^2) );

B_2 = 0;
for m = 1:l
    for n = m+1 : N
        B_2 = B_2 + 2*real(hat_c(n)*conj(v(n))* conj(hat_c(n-m))*v(n-m)); 
        %B_3 = B_3 + 2*real(conj(hat_c(n))*v(n)* hat_c(n-m)*conj(v(n-m))); 
        %B_2 = B_2 + 2*real(hat_c(n)*v(n)* conj(hat_c(n-m))*conj(v(n-m))); 
    end
end


B_est = B_1 + B_2;
end
