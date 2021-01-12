clc; clear;
seed = 7;
rng(seed);
Rep = 5000; % #MC replications
N = 125; 
k= 11;      %Number of regressors
beta = 0; %beta under H0 %VRAGEN


Z=normrnd(0,1, [N,k]);

alist = [0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0];
rholist = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 0.95];


Pz= Z*inv(Z.'*Z)*Z.';
Mz = eye(N)-Pz;

a = 1.5;
rho = 0.1;


%%%% QUESTION 2

rb_grid = 0:10:500;
LRbeta0 = zeros(length(rb_grid),Rep);
for rb=rb_grid
    rng(seed);
    for i=1:Rep
        
        %simulate psi's
        psi_1=chi2rnd(1);
        psi_k_1=chi2rnd(k-1);
        
        %Compute LR
        LRbeta0(find(rb_grid==rb), i)= 0.5*(psi_k_1 + psi_1-rb + sqrt((psi_k_1 + psi_1 + rb)^2 -4*rb*psi_k_1));  
    end
    
end
sortLR= sort(LRbeta0,2);
critLR= sortLR(:,round(0.95*Rep));

figure;
plot(rb_grid, critLR)
xlabel('rbeta_0')
ylabel('Critical values')
title('Critical value function of LR(beta) with k = 11')

%%%% QUESTION 1

rejectfreq2SLS = zeros(length(alist), length(rholist));
rejectfreqAR = zeros(length(alist), length(rholist));
rejectfreqscore = zeros(length(alist), length(rholist));
rejectfreqLR = zeros(length(alist), length(rholist));
% loop over a
for a = alist
    
    % loop over rho
    for rho = rholist
        rejectH02sls = 0;
        rejectH0AR = 0;
        rejectH0score = 0;
        rejectH0LR = 0;

        % 5000 MC simulations
        rng(seed);
        for i = 1:Rep
            % define model
            Pi = a * [1; zeros(10,1)];
            Sigma = [1, rho; rho , 1];

            Disturb  = mvnrnd([0, 0], Sigma, N);
            eps = Disturb(:,1);
            V = Disturb(:,2);
            X = Z*Pi + V;
            Y = X*beta + eps;
            
            
            
           
            % get 2SLS t statistic
            % stage 1
            Pihat = Z\X  ;

            % stage 2
            beta2sls = (Z*Pihat)\Y;

            s22sls= 1/(N-1)*(Y-X*beta2sls).'*(Y-X*beta2sls);
            varbeta2sls= s22sls*inv((X.'*Pz*X));

            t2sls= (beta2sls-beta)/sqrt(varbeta2sls);
            
            %AR
            AR = ((Y-X*beta).'*Pz*(Y-X*beta)/k)/((Y-X*beta).'*Mz*(Y-X*beta)/(N-k));
            
            % Score
            sighatee = 1/(N-k)*(Y-X*beta).'*Mz*(Y-X*beta);
            sighatev = 1/(N-k)*(Y-X*beta).'*Mz*X;
            rhohat = sighatev/sighatee.';
            Pitilde = inv(Z.'*Z)*Z.'*(X-eps*rhohat);
            score = 1/sighatee*(Y-X*beta).'*Z*Pitilde*inv(Pitilde.'*Z.'*Z*Pitilde)*Pitilde.'*Z.'*(Y-X*beta);
           
            % Compute the LR statistic
            kAR =(1/sighatee)*(Y-X*beta).'*Pz*(Y-X*beta);
            Sigvve = 1/(N-k)*X.'*Mz*X-(sighatev*sighatev)/sighatee;
            rbeta0= 1/Sigvve*Pitilde.'*Z.'*Z*Pitilde;
            LR=1/2*(kAR -rbeta0+sqrt((kAR+rbeta0)^2-4*rbeta0*(kAR-score)));
            
            % count rejections
            if abs(t2sls) > tinv(0.975,N-1);  
                rejectH02sls = rejectH02sls + 1;
            end
            if AR > chi2inv(0.95, k)/k 
                rejectH0AR = rejectH0AR +1;
            end
            if score > chi2inv(0.95,1) 
                rejectH0score = rejectH0score +1;
            end  
            if LR > critLR
                rejectH0LR = rejectH0LR +1;
            end
        end
        
        % fill matrices with amount of rejected
        rejectfreq2SLS(find(alist==a), find(rholist==rho)) = rejectH02sls/Rep;
        rejectfreqAR(find(alist==a), find(rholist==rho)) = rejectH0AR/Rep;
        rejectfreqscore(find(alist==a), find(rholist==rho)) = rejectH0score/Rep;
        rejectfreqLR(find(alist==a), find(rholist==rho)) = rejectH0LR/Rep;
    end

end


legendinfo=cell(length(a),1);
legendinfo{1}='a = 1';
for i=2:length(alist)
    legendinfo{i}=['a = ' num2str(alist(i))];
end

% plot figure for every rho 2 for stage least squares
figure;
plot(rholist, rejectfreq2SLS)
xlabel('rho')
ylabel('Rejection frequency')
title('Rejection frequency 2SLS t statistic')
legend(legendinfo)


% plot figure for every rho for AR test stat
figure;
plot(rholist, rejectfreqAR)
xlabel('rho')
ylabel('Rejection frequency')
title('Rejection frequency AR test statistic')
legend(legendinfo)


% plot figure for every rho for score test
figure;
plot(rholist, rejectfreqscore)
xlabel('rho')
ylabel('Rejection frequency')
title('Rejection frequency score test statistic')
legend(legendinfo)


% plot figure for every rho for LR test 
figure;
plot(rholist, rejectfreqLR)
xlabel('rho')
ylabel('Rejection frequency ')
title('Rejection frequency LR test statistic')
legend(legendinfo)


%%%% QUESTION 4
k_4 = 4;
rb_grid = 0:10:500;
LRbeta0 = zeros(length(rb_grid),Rep);
for rb=rb_grid
    rng(seed);
    for i=1:Rep
        %simulate psi's
        psi_1=chi2rnd(1);
        psi_k_1=chi2rnd(k_4-1);
        
        %Compute LR
        LRbeta0(find(rb_grid==rb), i)= 0.5*(psi_k_1 + psi_1-rb + sqrt((psi_k_1 + psi_1 + rb)^2 -4*rb*psi_k_1));  
    end
    
end
sortLR= sort(LRbeta0,2);
critLR= sortLR(:,round(0.95*Rep));

figure;
plot(rb_grid, critLR)
xlabel('rbeta_0')
ylabel('Critical values')
title('Critical value function of LR(beta) with k = 4')

%%%%%%%QUESTION 5%%%%%%%%%
%Question 5a

%clear ; clc;
load('assignmentweakinstruments.mat')

%define parameters
N = length(wage);
W = [ones(N,1) exper exper2 south smsa race];
Mw= eye(N)-W*inv(W.'*W)*W.';
Z = Mw* nearc2;
Y = Mw* wage;
X = Mw* ed;
Pz= Z*inv(Z.'*Z)*Z.';
Mz = eye(N)-Pz;


% get 2SLS t statistic
% stage 1
Pihat = Z\X  ;

% stage 2
beta2sls = (Z*Pihat)\Y;

s22sls= 1/(N-1)*(Y-X*beta2sls).'*(Y-X*beta2sls);
varbeta2sls= s22sls*inv((X.'*Pz*X));

Upper2sls = beta2sls + tinv(0.95, N-1)*sqrt(varbeta2sls);
Lower2sls = beta2sls - tinv(0.95, N-1)*sqrt(varbeta2sls);
fprintf( ' Upper bound: %.4f , lower bound: %.4f \n', Upper2sls, Lower2sls) 

%AR
AR = ((Y-X*beta).'*Pz*(Y-X*beta)/k)/((Y-X*beta).'*Mz*(Y-X*beta)/(N-k));

