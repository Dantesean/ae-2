clc; clear;
data = load('assignmentweakinstruments.mat');
W = [ones(3010,1) data.exper data.exper2 data.south data.smsa data.race];
MW = eye(size(W,1))-W*inv(W'*W)*W';
Y = MW*data.wage;
X = MW*data.ed;
PX = X*inv(X'*X)*X';
MX = eye(size(X,1))-PX;
Z = [data.nearc2 data.nearc4 data.nearc4a data.nearc4b];
Z = MW*Z;
PZ = Z*inv(Z'*Z)*Z';
MZ = eye(size(Z,1))-PZ;
N = size(Y,1);
k = size(Z,2);
beta_0=(-5:0.01:5);
psi_1 = chi2rnd(1);
psi_k = chi2rnd(k-1);


    for j=1:length(beta_0) 
                b_2sls = inv(X'*PZ*X)*X'*PZ*Y;
                res_2sls = Y-X*b_2sls;
                var_2sls = res_2sls'*res_2sls/(N-1)*inv(X'*PZ*X);
                se_2sls = sqrt(var_2sls);
                
                error_0 = Y-X*beta_0(j);
                sigma_hat_ee = (1/(N-k))*error_0'*MZ*error_0;
                sigma_hat_eV = (1/(N-k))*error_0'*MZ*X;
                sigma_hat_VV = (1/(N-k))*X'*MZ*X;
        
                rho_hat = sigma_hat_eV/sigma_hat_ee;
                pi_tilde = inv(Z'*Z)*Z'*(X-error_0*rho_hat);
                Zpi_tilde = Z*pi_tilde;
                PZpi = Zpi_tilde*inv(Zpi_tilde'*Zpi_tilde)*Zpi_tilde';

                sigma_hat_VVe = sigma_hat_VV - (sigma_hat_eV^(2)/sigma_hat_ee);
        
                r_beta(j) = (1/sigma_hat_VVe)*Zpi_tilde'*Zpi_tilde;
                %psi_1 = chi2rnd(1);
                %psi_k = chi2rnd(k-1);
                
                t_2sls(j) = (b_2sls-beta_0(j))/sigma_hat_ee;
                AR(j) = ((error_0'*PZ*error_0)/k)/sigma_hat_ee;
                LM(j) = ((error_0'*PZpi*error_0))/sigma_hat_ee;
                LR(j) = 0.5*(k*AR(j) - r_beta(j) + sqrt((k*AR(j)+r_beta(j))^2 - 4*r_beta(j)*(k*AR(j)-LM(j))));
                LR_crit(j) = 0.5*(psi_k + psi_1 - r_beta(j) + sqrt((psi_k + psi_1 + r_beta(j))^2-4*r_beta(j)*psi_k));
    end
    


%%%% 5e
figure
plot(beta_0,t_2sls, 'LineWidth', 1.5)
hold on
plot(beta_0,AR, 'LineWidth', 1.5)
hold on
plot(beta_0,repmat(chi2inv(0.95,k)/k,length(beta_0)), 'LineWidth', 1.5)
hold on
plot(beta_0, repmat(norminv(0.95), length(beta_0)), 'LineWidth', 1.5)

figure
plot(beta_0,LM, 'LineWidth', 1.5)
hold on
plot(beta_0,LR, 'LineWidth', 1.5)
hold on
plot(beta_0,repmat(chi2inv(0.95,1),length(beta_0)), 'LineWidth', 1.5)
hold on
plot(beta_0,LR_crit, 'LineWidth', 1.5)
