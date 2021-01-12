clc;
clear;

%%%% Start time
tic

% Parameter values
beta = 0;                                  % parameter of interest (single regressor)
rho  = [(0:0.1:0.7) 0.9 0.95];             % degree of endogeneity
k = 11;                                    % number of instruments
a =[1.5 0.7 0.5 0.3 0.15 0.07 0.03 0];     % instrument strength
N = 125;                                   % sample size
e11=[1;zeros(10,1)];

% Exogenous instruments, (keep) fixed in repeated samples
Z = normrnd(0,1,N,k);

% Creating matrices to collect the results of Reps replications
Reps = 5000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_beta=0:0.01:500;
LR_crit=LRcrit(r_beta,11,Reps);

figure(5);
plot(r_beta,LR_crit)
xlabel('r(beta)')
ylabel('95% critical value')
ylim([0 20])
title('Critical values of LR-statistic for k=11')


%%%% End time
toc