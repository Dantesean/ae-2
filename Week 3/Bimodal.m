%  Bimodal.m
%
%  This programme can be used to generate artificial data that are
%  bimodally distributed and then estimated nonparametrically
%
% (c) K.J. van Garderen
%
%  You can use F9 to run a selection of the program to run it in steps
%
% Create Mixed Normal/Gaussian bivariate distribution Object
%   using   mixG_obj = gmdistribution(mu,sigma,mixp);
clear;
mu    = [1 1;0 -2];                       % mean  of first and second bivariate Gaussian
sigma = cat(3,[2 0;0 .5],[1 0;0 1]);      % covariance of first and second bivariate Gaussian
mixp  = [0.9,0.9];                        % mixing proportions (with mixp = [0.1,0.9]; the second mode is much more important

mixG_obj = gmdistribution(mu,sigma,mixp); % This is the mixed Gaussian bivariate distribution with 2 modes


%  set the number of observations nobs
nobs = 1000;

% draw sample of nobs(ervations)  from this bivariate distribution
X = random(mixG_obj,nobs);
                                % matlab is case sensitive so X is different from x
                                % X(:,1) takes the first column in X  
                                % X(2,:) takes the 2nd row of X

% Define the midpoints at which we estimate the density
 xmid = linspace(-6,6,100)';     % 100 linearly space values between -6 and 6                              

% Estimate the density using your own kernel density estimator 
% using midpoints xmid and let the program determine the bandwidth 
  [Xpoints,pdf_est1,bndw] = npdensity_kjvg(X(:,1),xmid,0);
 sprintf('Bandwidth used = %6.3g',10);

  % graph the result
  plot(Xpoints,pdf_est1);
  hold on;            % plot everything in 1 graph
  [Xpoints,pdf_est2,bndw] = npdensity_kjvg(X(:,2),xmid,0);
  plot(Xpoints,pdf_est2);
  
% check that the estimated density integrates to (almost) 1 over the xmid
% range
  binsize = (Xpoints(2)-Xpoints(1)); % since equally spaced 
  sprintf('The density integrates to %6.3f with bandwidth %6.3g and binsize %6.3g ',sum(pdf_est1).*binsize,bndw,binsize)

  
% use matlabs own kerneldensity routine to compare with your result
k_x1 = ksdensity(X(:,1),Xpoints);   % this is MATLAB's own kernel density routine
plot(Xpoints,[pdf_est1, k_x1]);
hold on;
k_x2 =ksdensity(X(:,2),Xpoints);  
plot(Xpoints,[pdf_est2, k_x2]);

% you can plot the true bivariate distribution
% use it to explain why one marginal distribution looks much more bimodal
% than the other
ezsurf(@(x,y)pdf(mixG_obj,[x y]),[-6 6],[-6 6])


