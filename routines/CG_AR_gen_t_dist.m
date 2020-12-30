function c_norm=CG_AR_gen_t_dist(N,p,sigma2_c,lambda,scale)

p_abs = abs(p);

Ntrans=100*ceil(3*log(10)/abs( log( max(p_abs) ) ) ); 

K = N+Ntrans;

% Innovations
innov = sqrt(1/2)*(randn(1,K)+1j.*randn(1,K));

rho =  poly(p);

cll = filter(1,rho,innov);
R = gamrnd(lambda,scale,1,1);
cl = sqrt(1./R).*cll;

c = cl(1+Ntrans:end).';

sigma2_c_est = mean(abs(c).^2);

c_norm = (c./sqrt(sigma2_c_est))*sqrt(sigma2_c);

end