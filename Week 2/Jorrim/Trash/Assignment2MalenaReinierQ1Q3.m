clc;
clear;

%%%% Start time
tic

% Parameter values
beta = 0;                                  % parameter of interest (single regressor)
rho  = [(0:0.1:0.7) 0.9 0.95];                 % degree of endogeneity
k = 11;                                    % number of instruments
a =[1.5 0.7 0.5 0.3 0.15 0.07 0.03 0];     % instrument strength
N = 125;                                   % sample size
e11=[1;zeros(10,1)];

% Exogenous instruments, (keep) fixed in repeated samples
Z = normrnd(0,1,N,k);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1 & 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Used to store the rejection frequencies
reject_mat = zeros(length(a),length(rho)); 
reject_mat_AR = zeros(length(a),length(rho)); 
reject_mat_LM = zeros(length(a),length(rho)); 
reject_mat_LR = zeros(length(a),length(rho)); 


% Creating matrices to collect the results of Reps replications
Reps = 5000;


 for i=1:length(a)
     pi = a(i)*e11;
     
     for j=1:length(rho)
        
        rng(11045301);
     
        reject = 0;
        reject_AR = 0;
        reject_LM = 0;
        reject_LR = 0;
        
         for s=1:Reps
            errors= mvnrnd([0;0],[1 rho(j); rho(j) 1],N);
            eps=errors(:,1);
            v=errors(:,2);
            x = Z*pi + v;
            y = x*beta + eps;
            
            pi_hat=(Z\x);
            xhat = Z*pi_hat;
            b_2SLS =  xhat\y;      % calculated as  TSLS
            e = y - x*b_2SLS; 
            t_2SLS = (b_2SLS - beta)/sqrt(e'*e/(N*(xhat'*xhat))); 
            
            if abs(t_2SLS)>tinv(0.975,N-1)
                reject=reject+1;
            end
            
            % Score LM statistic (question 3) (slides 54-59)
            s_eps = 1/(N-k)*y'*(eye(N)-Z*inv(Z'*Z)*Z')*y;
            s_epsV = 1/(N-k)*y'*(eye(N)-Z*inv(Z'*Z)*Z')*x;
            rho_hat = s_epsV/s_eps;
            ZPi = Z*inv(Z'*Z)*Z'*(x-y*rho_hat);
            LM_stat = 1/s_eps*y'*(ZPi*inv(ZPi'*ZPi)*ZPi')*y;
            
            if LM_stat>chi2inv(0.95,1)
                reject_LM=reject_LM+1;
            end
            
            % AR statistic (question 3) (slide 51)
            
            %AR_stat = (N-k)/k*(y'*Z*inv(Z'*Z)*Z'*y)/(y'*(eye(N)-Z*inv(Z'*Z)*Z')*y);
              
            P_z=Z*inv(Z'*Z)*Z';
            M_z=eye(N)-P_z;
            yPy=(eps'*P_z*eps)/k;
            yMy=(eps'*M_z*eps)/(N-k);
            AR_stat= yPy/yMy;
            
            if AR_stat>chi2inv(0.95,k)/k
            reject_AR=reject_AR+1;
            end
            
            % LR statistic (question 3) (slides 63-66)
            s_VV = 1/(N-k)*x'*(eye(N)-Z*inv(Z'*Z)*Z')*x;
            s_VVeps = s_VV - s_epsV^2/s_eps;
            r_beta0 = 1/s_VVeps*(inv(Z'*Z)*Z'*(x-y*rho_hat))'*Z'*Z*inv(Z'*Z)*Z'*(x-y*rho_hat);
            LR_stat = 1/2*(k*AR_stat - r_beta0 + sqrt((k*AR_stat + r_beta0)^2 - 4*r_beta0*(k*AR_stat - LM_stat)));
            
            if LR_stat>LRcrit(r_beta0,k,Reps)
                reject_LR=reject_LR+1;
            end  
            
         end
       
         reject_mat(i,j)=reject/Reps;
         reject_mat_AR(i,j)=reject_AR/Reps;
         reject_mat_LM(i,j)=reject_LM/Reps;
         reject_mat_LR(i,j)=reject_LR/Reps;
         
     end
  
 end

figure(1);
for i=1:length(a)
    plot(rho, reject_mat(i,:))
    hold on
end

legendinfo=cell(length(a),1);
legendinfo{1}='a = 1.5';
legendinfo{2}='a = 0.7';
legendinfo{3}='a = 0.5';
legendinfo{4}='a = 0.3';
legendinfo{5}='a = 0.15';
legendinfo{6}='a = 0.07';
legendinfo{7}='a = 0.03';
legendinfo{8}='a = 0';

legend(legendinfo,'location','northwest')
title(sprintf('T-statistic'))
xlabel('Rho')
ylabel('Rejection frequency')
hold off

figure(2);
for i=1:length(a)
    plot(rho, reject_mat_LM(i,:))
    hold on
end

legendinfo=cell(length(a),1);
legendinfo{1}='a = 1.5';
legendinfo{2}='a = 0.7';
legendinfo{3}='a = 0.5';
legendinfo{4}='a = 0.3';
legendinfo{5}='a = 0.15';
legendinfo{6}='a = 0.07';
legendinfo{7}='a = 0.03';
legendinfo{8}='a = 0';

legend(legendinfo,'location','northwest')
title(sprintf('LM-statistic'))
xlabel('Rho')
ylabel('Rejection frequency')
hold off

figure(3);
for i=1:length(a)
    plot(rho, reject_mat_AR(i,:))
    hold on
end

legendinfo=cell(length(a),1);
legendinfo{1}='a = 1.5';
legendinfo{2}='a = 0.7';
legendinfo{3}='a = 0.5';
legendinfo{4}='a = 0.3';
legendinfo{5}='a = 0.15';
legendinfo{6}='a = 0.07';
legendinfo{7}='a = 0.03';
legendinfo{8}='a = 0';

legend(legendinfo,'location','northwest')
title(sprintf('AR-statistic'))
xlabel('Rho')
ylabel('Rejection frequency')
hold off

figure(4);
for i=1:length(a)
    plot(rho, reject_mat_LR(i,:))
    hold on
end

legendinfo=cell(length(a),1);
legendinfo{1}='a = 1.5';
legendinfo{2}='a = 0.7';
legendinfo{3}='a = 0.5';
legendinfo{4}='a = 0.3';
legendinfo{5}='a = 0.15';
legendinfo{6}='a = 0.07';
legendinfo{7}='a = 0.03';
legendinfo{8}='a = 0';

legend(legendinfo,'location','northwest')
title(sprintf('LR-statistic'))
xlabel('Rho')
ylabel('Rejection frequency')
hold off


%%%% End time
toc
