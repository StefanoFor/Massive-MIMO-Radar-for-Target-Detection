clear;
close all
clc

%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Monte Carlo runs
MC=10^5;


% Definition of the complex coefficients of the AR(p) disturbance process
nu1 = -0.4;
nu2 = -0.2;
nu3 = 0;
nu4 =  0.1;
nu5 = 0.3;
nu6 = 0.35;
p(1) = 0.5*exp(1j*2*pi*nu1);
p(2) = 0.6*exp(1j*2*pi*nu2);
p(3) = 0.7*exp(1j*2*pi*nu3);
p(4) = 0.4*exp(1j*2*pi*nu4);
p(5) = 0.5*exp(1j*2*pi*nu5);
p(6) = 0.6*exp(1j*2*pi*nu6);
rho = - poly(p);

% Power of the disturbance process
sigma2_c = 1;

% Parameters of the t-distributed innovations
lambda = 2;
eta = lambda/(lambda-1);
scale=eta/lambda;

% Number of antennas
Nvect= floor(logspace(2,4,20));    % Number of antennas
Nl=length(Nvect);

% Nominal PFA and related threshold
PFA_nom = 10^(-4);
th = chi2inv(1-PFA_nom,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial parameter of the Cell Under Test (CUT)
CUT = -0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR = -20;
sigma2_b=(sigma2_c)*10^(SNR/10);   % Power of complex amplitude
alpha = sqrt(sigma2_b)*exp(1j*2*pi*rand(1));  % complex target amplitude

PD_W_T_est = zeros(1,Nl);
PD_nom = zeros(1,Nl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for il=1:Nl
    tic;
    PD_W_T=0;
    
    N = Nvect(il)
    l=floor(N^(1/4));
    m=[0:N-1].';
    v=exp(1j*2*pi*CUT*m);   % Steering vector
    
    n2 = (v'*v);
    
    parfor ins=1:MC
        
        x_H0 = AR_gen_t_dist(N,p,sigma2_c,lambda,scale);
        x = alpha * v + x_H0;
        
        % Estimation of alpha under H1
        hat_alpha = v'*x/n2;
        
        % Wald Test
        W_T = Wald_test(N,v,x,hat_alpha,l,n2);
        
        % PD Wald test
        PD_W_T=PD_W_T+(sign(W_T-th)+1)/2;
        
        B_est_mod_0(ins) = B_matrix_est(N,v,x_H0,0,l);
        
    end
    
    toc/60
    
    B_H0_mean = mean(B_est_mod_0);
    
    PD_W_T_est(il)=PD_W_T/MC;
    
    non_cen_par = 2*abs(alpha)^2*n2^2/B_H0_mean;
    PD_nom(il) = marcumq(sqrt(non_cen_par),sqrt(th));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save PD_CUT_meno02_AR6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

figure;
semilogx(Nvect,PD_W_T_est,'LineWidth',2)
hold on
semilogx(Nvect,PD_nom,'LineWidth',2)
%axis([10 1000 10^(-2) 10^(-1)])
grid on;
legend('Wald Test','Nominal');
xlabel('N');
ylabel('Probability of Detection (P_{D})')


