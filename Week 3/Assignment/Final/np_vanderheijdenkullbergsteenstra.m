%%
clc; clear all; clf;

% make data
TempCollision = importdata('collisions_grouped_3_no_rain.csv').data; % only load in no rain data
temp = TempCollision(:,1); % average temperature over all collisions on that day
collision = TempCollision(:,2); % the amount of collisions on that day

% basic data representation
% f1 = figure(1);
% ax1 = axes('Parent', f1);
% histogram(ax1, temp,50, 'FaceColor', 'b');
% title(ax1, 'Histogram of temperatures at collisions without rain');
% xlabel(ax1, 'Temperature');
% ylabel(ax1, 'Frequency');
% 
% f2 = figure(2);
% ax2 = axes('Parent', f2);
% histogram(ax2, collision,50, 'FaceColor', 'b');
% title(ax2, 'Histogram of collisions');
% xlabel(ax2, 'Amount of collisions');
% ylabel(ax2, 'Frequency');
% 
% f3 = figure(3);
% ax3 = axes('Parent', f3);
% yyaxis left
% plot(ax3, temp, 'b')
% ylabel(ax3, 'Temperature')
% yyaxis right
% plot(ax3, collision, 'r')
% title(ax3, 'Temperatures/Collisions over time')
% ylabel(ax3, 'Amount of collissions')
% xlabel(ax3, 'Time');
% ax3.YAxis(1).Color = 'b';
% ax3.YAxis(2).Color = 'r';

%% OLS

X_ols = [ones(length(temp),1) temp temp.^2 temp.^3 temp.^4]; % add intercept and powers
[b,bint,res,res_int,stats] = regress(collision,X_ols); % ols
% some other statistic, could be used for analysis
R_sq = stats(1);
F_stat = stats(2);
F_pvalue = stats(3);
s2_hat = stats(4);


%% Kernel density
temp_midpoints = linspace(min(temp),max(temp),100)'; % create midpoints

% k_pdf is the calculated density, according to the method specified
[Xmidpoints_used,k_pdf,bandwidth_used] = npdensity_DvdHWKWS(temp, temp_midpoints, -1); % change -2 to what you want
% -1: Epanechnikov, -2: Silverman, >0: manual entry, else Gaussian
[Xmidpoints_used_silverman,k_pdf_silverman,bandwidth_used_silverman] = npdensity_DvdHWKWS(temp, temp_midpoints, -2); %-2 is silverman
[Xmidpoints_used_gaus,k_pdf_gaus,bandwidth_used_gaus] = npdensity_DvdHWKWS(temp, temp_midpoints, 0); %-2 is silverman

binsize = (temp_midpoints(2)-temp_midpoints(1));
eq1 = sum(k_pdf).*binsize; % check if it equals 1

% matlab's function
k_pdf_matlab=ksdensity(temp,temp_midpoints);


%% Kernel regression
% first output is the x-midpoints used, 
% second output is the estimation of the non-parametric regression, 
% third output is used bandwidth
% fourth output is sigma^2 hat estimation
% fifth and sixth output are the confidence intervals

% first call npregress function
[X_mid_gaus,m_gaus,bandwidth_gaus,s_2_gaus, CI_low_gaus, CI_high_gaus] = npregress_DvdHWKWS(collision,temp,temp_midpoints,0); % use gaussian
% specify method: -1: Epanechnikov, -2: Silverman, -3: CV, >0: manual entry, else Gaussian


% % same for CV
[X_mid_CV,m_CV,bandwidth_CV,s_2_CV, CI_low_CV, CI_high_CV] = npregress_DvdHWKWS(collision,temp,temp_midpoints,-3); % now cv


[X_mid_silverman,m_silverman,bandwidth_silverman,s_2_silverman, CI_low_silverman, CI_high_silverman] = npregress_DvdHWKWS(collision,temp,temp_midpoints,-2); % now cv



%% calculate ols prediction
X_mid_ols = [ones(length(X_mid_gaus),1) X_mid_gaus X_mid_gaus.^2 X_mid_gaus.^3 X_mid_gaus.^4];
y_hat = X_mid_ols*b;

y_hat_se = s2_hat*(X_mid_ols*inv(X_mid_ols'*X_mid_ols)*X_mid_ols' + eye(length(X_mid_ols)));
CI_ols_low = y_hat-1.96*sqrt(diag(y_hat_se));
CI_ols_high = y_hat+1.96*sqrt(diag(y_hat_se));

%% plots
% Plot 1 & 5
f1 = figure(1);
ax1 = axes('Parent', f1);
histogram(ax1, (temp/max(temp)*100),50, 'FaceColor', 'b', 'Normalization', 'probability');
hold on 

plot(ax1, k_pdf_matlab, 'c')
title('PDF estimation')
hold on

plot(ax1, k_pdf_silverman, 'r');
hold on

plot(ax1, k_pdf_gaus, 'b');
legend(ax1, 'True data', 'Matlabs own density estimator', 'Silverman', 'Gaussian', 'Location','NorthEast')
title(ax1, 'Histogram and PDF estimation')

% Plot 4 & 6
f2 = figure(2);
ax2 = axes('Parent', f2);
scatter(ax2, temp,collision, 'y');
title(ax2, 'Scatterplot and OLS 4th power')
xlabel(ax2, 'Temperature')
ylabel(ax2, 'Amount of collisions')
hold on
plot(ax2, X_mid_gaus,y_hat, 'c')
hold on
plot(ax2, X_mid_gaus, CI_ols_low, 'c--')
plot(ax2, X_mid_gaus, CI_ols_high, 'c--', 'HandleVisibility', 'off')
hold on

plot(ax2, X_mid_gaus,m_gaus, 'k-', 'HandleVisibility', 'off')
hold on
plot(ax2, X_mid_gaus, CI_low_gaus, 'k--', 'HandleVisibility', 'off')
plot(ax2, X_mid_gaus, CI_high_gaus, 'k--', 'HandleVisibility', 'off')
hold on

% same for CV
plot(ax2, X_mid_CV,m_CV, 'r-')
plot(ax2, X_mid_CV, CI_low_CV, 'r--')
plot(ax2, X_mid_CV, CI_high_CV, 'r--', 'HandleVisibility', 'off')
hold on

plot(ax2, X_mid_silverman,m_silverman, 'b-')
hold on
plot(ax2, X_mid_silverman, CI_low_silverman, 'b--')
plot(ax2, X_mid_silverman, CI_high_silverman, 'b--', 'HandleVisibility', 'off')

% we dont use epanechnikov
% plot(ax2, X_mid_epanechnikov,m_epanechnikov, 'g-')
% hold on
% plot(ax2, X_mid_epanechnikov, CI_low_epanechnikov, 'c--')
% plot(ax2, X_mid_epanechnikov, CI_high_epanechnikov, 'c--', 'HandleVisibility', 'off')

title(ax2, 'Regression estimates')

legend(ax2, 'Scatterplot', 'OLS', 'OLS CI', 'CV', 'CV CI', 'Gaussian/Silverman', 'Gaussian/Silverman CI', 'Location', 'NorthEast')

%% functions

function [Xmidpoints_used,k_pdf,bandwidth_used] = npdensity_DvdHWKWS(X,Xmidpoints,bandwidth)

% INPUTS
    %  X          : Nx1  vector of data  (if X contains more columns only the
    %                                      first column will be used)
    %  Xmidpoints : nrbins x 1 columnvector of bin midpoints
    %  bandwidth  :  if > 0 then this is the bandwidth that is used in the kernel density estimation   
    %                (if -1 , let the program determine the bandwidth based
    %                on Epanechnikov), but we don't use this one
    %                if -2, based on Silverman
    %                for CV density we also need estimation of regression,
    %                so that's implemented in the regress function
    %                else, based on Gaussian
% OUTPUTS
    % Xmidpoints_used     : the midpoints might be defined in this function if the
    %                  input = 0
    % k_pdf          : Density evaluated at Xmidpoints_used
    % bandwidth_used : bandwidth used
% local variables used
    % firstb lastb binsize nrbins bandwidth_used Nbandwidth_used k_pdf DD 
    % I J nr makebins mkbandwidth  meanX stdX Xb  Z KX;

   [nr,nc]  = size(X);
   
   X = X(:,1);         % take only first column of X
   meanX = mean(X);
   stdX  = std(X)';
   iota = ones(nr,1); 
    X = sortrows(X);   % order the observations for column 1
   
    
  %  nrbins is the number of bins/midpoints/gridpoints where the density is calculated  
    if Xmidpoints==0
        nrbins = 20;         % our DEFAULT number of bins (evaluation points) when it is NOT user defined
        Xmidpoints_used = linspace(X(floor(0.01*nr+1),1),X(floor(0.99*nr),1),nrbins)';  
        % this creates an equally spaced set of X midpoints between the 1 and 99 percentile 
    else
        [nrbins,ncbins] = size(Xmidpoints);
        Xmidpoints_used = Xmidpoints;
    end
    
    firstb = Xmidpoints_used(1,:);
    lastb  = Xmidpoints_used(nrbins,:);
    binsize = ((lastb - firstb)/nrbins); % average binsize in this case

    delta =0.7764;               % see (9.11) and table 9.1 Cameron & Trivedi
    if bandwidth>0 % use the bandwidth given in the proc argument
        bandwidth_used = bandwidth;
    elseif bandwidth == -1 % Epanechnikov, we don't use Epanechnikov
        bandwidth_used = 2.345*stdX*nr^(-1/5);
    elseif bandwidth == -2 % Silverman
        bandwidth_used = 1.3643*delta*nr^(-1/5)*min(stdX,iqr(X)/1.349);
    else                   % use Gaussian
        bandwidth_used = 1.3643*delta*stdX*nr^(-1/5); % Gaussian


    end    % endif 

k_pdf  = zeros(nrbins,1);  %  this will contain the density estimates at the midpoints_used
                           %  bin by bin calculations since reduces the amount of workspace used vs k_pdf at each observation %
  
    for J=1:nrbins         % for each bin       %  
       Xb = Xmidpoints_used(J,1);          % one bin at a time  %
       Z = (iota*Xb - X)/bandwidth_used;
       if bandwidth == -1
           KX = max([3/4*(1-Z.^2), zeros(nr,1)]') / bandwidth_used; % Epanechnikov Kernel
       else
           KX = pdf('Normal',Z,0,1)/bandwidth_used;     % rest uses Gaussian Kernel%
       end
       k_pdf(J,1) = mean(KX);
    end    % for; 

end


function [Xmidpoints_used,m_regress,bandwidth_used,s_2, CI_low, CI_high] = npregress_DvdHWKWS(Y,X,Xmidpoints,bandwidth)

% INPUTS
    %  X          : Nx1  vector of data  (if X contains more columns only the
    %                                      first column will be used)
    %  Xmidpoints : nrbins x 1 columnvector of bins. If ==0 then use default
        %  bandwidth  :  if > 0 then this is the bandwidth that is used in the kernel density estimation   
    %                (if -1 , let the program determine the bandwidth based
    %                on Epanechnikov), but we don't use this one
    %                if -2, based on Silverman
    %                if -3, use CV density
    %                else, based on Gaussian
% OUTPUTS
    % Xmidpoints_used     : the midpoints might be defined in this function if the
    %                  input = 0
    % k_pdf          : Density evaluated at Xmidpoints_used
    % bandwidth_used : bandwidth used
% local variables
    % firstb lastb binsize nrbins bandwidth_used Nbandwidth_used k_pdf DD 
    % I J nr makebins mkbandwidth  meanX stdX Xb  Z KX;

   [nr,nc]  = size(X);
   YX = horzcat(Y,X);    % keep Xi and Yi together when sorting
   YX = sortrows(YX,2);  %  sort observations according to the X variable (second column)
   
   Y = YX(:,1);
   X = YX(:,2);   % take only first column of X to determine the largest and smallest values
  
   meanX = mean(X);
   stdX  = std(X)';
   iota  = ones(nr,1); 
  %  X = sortrows(X);   % order the observations for column 1
   
    
  %  nrbins is the number of bins/midpoints/gridpoints where the density is calculated  
    if Xmidpoints==0       % then use default
        nrbins = 100;         %our default number of bins (evaluation points) when it is not user defined
        Xmidpoints_used = linspace(X(floor(0.01*nr+1),1),X(floor(0.99*nr),1),nrbins)';  
        % this creates an equally spaced set of X midpoints between the
        %                                         1 and 99 percentile  of X
    else
        [nrbins,ncbins] = size(Xmidpoints);
        Xmidpoints_used = Xmidpoints;
    end
    
    firstb = Xmidpoints_used(1,:);
    lastb  = Xmidpoints_used(nrbins,:);
    binsize = ((lastb - firstb)/nrbins); % average binsize in this case

    delta =0.7764;               % see (9.11) and table 9.1 Cameron & Trivedi
    if bandwidth>0 % use the bandwidth given in the proc argument
        bandwidth_used = bandwidth;
    elseif bandwidth == -1 % Epanechnikov don't use Epanechnikov
        bandwidth_used = 2.345*stdX*nr^(-1/5);
    elseif bandwidth == -2 % Silverman
        bandwidth_used = 1.3643*delta*nr^(-1/5)*min(stdX,iqr(X)/1.349);
    elseif bandwidth == -3 % CV
        h = linspace(0.2,2,41);
        errors =  zeros(length(X),1);
        CV = zeros(length(h),1);
        for i = 1:length(h)
            h_val = h(i);
            for j = 1:length(X)
                xj = X(j);
                yj = Y(j);
                Z = (X - iota*xj)/h_val;
                KX = pdf('Normal',Z,0,1);
                wiih = pdf('Normal',0,0,1)/(sum(KX));
                YKX = Y.*KX;
                m_reg = (sum(YKX)-yj*wiih)/(sum(KX)-wiih);

                errors(j) = ((yj-m_reg)/(1-wiih))^2;
            end
            CV(i) = mean(errors.^2);
        end
        bandwidth_used = h(find(CV == min(CV))); % min value of CV
    else                   % use Gaussian
        bandwidth_used = 1.3643*delta*stdX*nr^(-1/5); % Gaussian

    end    %endif 


    K_f = zeros(nrbins,1);
    m_regress= zeros(nrbins,1); % for each bin
    s_2 = zeros(nrbins,1);
      %  bin by bin due to workspace limitations %
  
    for J=1:nrbins                       % for each bin       %
       Xb = Xmidpoints_used(J,1);          % one bin at a time  %
       Z = (iota*Xb - X)/bandwidth_used;
       if bandwidth == -1
           KX = max([3/4*(1-Z.^2), zeros(nr,1)]')' / bandwidth_used; % Epanechnikov Kernel
       else
           KX = pdf('Normal',Z,0,1)/bandwidth_used;     % rest uses Gaussian Kernel%
       end
      
       YKX= Y.*KX;
       m_regress(J,1) = mean(YKX)/mean(KX);
    end    % for; 
    
    % calculate weights
    pi_i = zeros(nr,nrbins);
    for J = 1:nrbins
        pi_i(:,J) = pdf('Normal',(X - Xmidpoints_used(J))/bandwidth_used,0,1); % calc weights
    end
    pi_r = pi_i./sum(pi_i);  % relative weights (make them sum up to 1)

    % calculate errors
    eps = zeros(nr,1);
    for j = 1:nr
            xj = X(j);
            yj = Y(j);
            Z = (X - iota*xj)/bandwidth_used;
            KX = pdf('Normal',Z,0,1);
            YKX = Y.*KX;
            m_reg = sum(YKX)/sum(KX);
            eps(j) = yj - m_reg;
    end
    
    % calculate estimator of variance
    for J = 1:nrbins
        s_2(J) = pi_r(:,J)'*eps.^2;
    end
    
    %last term for CI
    for J=1:nrbins                       % for each bin       %
       Xb = Xmidpoints_used(J,1);          % one bin at a time  %
       Z = (X - iota*Xb)/bandwidth_used;
       if bandwidth == -1
           K_sq = @(z)(3/4*(1-z.^2)).^2;
           int_K = integral(K_sq,-1,1);
           f_0 = mean(max([3/4*(1-Z.^2), zeros(nrows,1)]')/bandwidth_used)/bandwidth_used; % Epanechnikov Kernel
       else
           K_sq = @(z)(1/sqrt(2*pi)*exp(-1/2*z.^2)).^2;
           int_K = integral(K_sq,-Inf,Inf);
           f_0 = mean(pdf('Normal',Z,0,1)/bandwidth_used);
       end
       
    K_f(J) = sqrt(int_K / (f_0*nr*bandwidth_used));
    end    % for; 

    
    % bootstrap bias, code from week 1
    bootrep = 100;
    boot_avg = zeros(nrbins,1);
    for J = 1:nrbins
        m_boot = zeros(bootrep,1);
        for r=1:bootrep
            index=unidrnd(nr,nr,1);    % randomly selected
        
            XB=X(index);              
            YB = Y(index);
        
            Xb = Xmidpoints_used(J,1);
            Z = (iota*Xb - XB)/bandwidth_used;
            KX = pdf('Normal',Z,0,1)/bandwidth_used;     
        
            YKX= YB.*KX;
            m_boot(r,1) = mean(YKX)/mean(KX);
        end
        boot_avg(J,1) = mean(m_boot);
    end
    
    bias = m_regress - boot_avg;
    
    % calculate naive/simple CI for regression
    % CI_low = m_reg-1.96*sqrt(s_2);
    % CI_high = m_reg+1.96*sqrt(s_2);
    
    
    % bias-corrected CI estimate would be
    CI_low = m_regress-bias-1.96*sqrt(s_2).*K_f;
    CI_high = m_regress-bias+1.96*sqrt(s_2).*K_f;
    
end

