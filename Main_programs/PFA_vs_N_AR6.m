
clear all;
close all
clc

%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Monte Carlo runs
MC=10^6;

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

% Normalized PSD of the AR(p) disturbance
plot_normalized_PSD(rho);

% Power of the disturbance process
sigma2_c = 1;

% Parameters of the t-distributed innovations
lambda = 2;
eta = lambda/(lambda-1);
scale=eta/lambda;

% Number of antennas
Nvect= floor(logspace(2,4,20));    % Number of antennas
Nvect= [Nvect,21544,46415,100000];

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

PFA_W_T_est = zeros(1,Nl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for il=1:Nl
    tic;
    
    PFA_W_T=0;
    
    N = Nvect(il)
    l=floor(N^(1/4));
    m=[0:N-1].';
    v=exp(1j*2*pi*CUT*m);   % Steering vector
    
    n2 = (v'*v);
    
    parfor ins=1:MC
        
        x = AR_gen_t_dist(N,p,sigma2_c,lambda,scale);
        
        % Estimation of the signal parameter alpha
        hat_alpha = v'*x/n2;
        
        % Wald Test
        W_T = Wald_test(N,v,x,hat_alpha,l,n2);
        
        % PFA Wald test
        PFA_W_T=PFA_W_T+(sign(W_T-th)+1)/2;
        
    end
    
    toc/60
    
    PFA_W_T_est(il)=PFA_W_T/MC;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save PFA_CUT_meno02_AR6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

time_total = toc/60

figure;
loglog(Nvect,PFA_W_T_est,'LineWidth',2)
hold on
loglog(Nvect,PFA_nom*ones(1,Nl),'LineWidth',2)
%axis([Nvect(1) Nvect(end) 10^(-5) 10^(-2)])
grid on;
legend('Wald test','Nominal');
xlabel('N');
ylabel('Probability of False Alarm (P_{FA})')


