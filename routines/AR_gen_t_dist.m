function c_norm=AR_gen_t_dist(N,p,sigma2_c,lambda,scale)

p_abs = abs(p);

Ntrans=100*ceil(3*log(10)/abs( log( max(p_abs) ) ) ); 

K = N+Ntrans;

% Innovations
w = sqrt(1/2)*(randn(1,K)+1j.*randn(1,K));
R = gamrnd(lambda,scale,1,K);
innov = sqrt(1./R).*w;

rho =  poly(p);

cl = filter(1,rho,innov);
c = cl(1+Ntrans:end).';

sigma2_c_est = mean(abs(c).^2);

c_norm = (c./sqrt(sigma2_c_est))*sqrt(sigma2_c);

end



