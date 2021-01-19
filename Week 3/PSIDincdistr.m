%   PSIDincdistr.m
%
% This programme uses the Panel Study of Income Dynamics to estimate the income disitribution non-parametrically
% It highlights some of the practical issues with real data and you're invited to think about solutions
%

clear;
PSID =importdata('indearnPSID.asc');
earn = PSID.data;

earn1mid = linspace(0,200000,100)';     % 100 linearly space values between 0 and 200K  
[earnpoints,pdf_est,bndw] = npdensity_kjvg(earn,earn1mid,0);

plot(earnpoints,pdf_est);


% Check that the approximate density integrates to 1 and check if we miss
% important part of the data
binsize = (earnpoints(2)-earnpoints(1)); % since equally spaced 
sprintf('The density integrates to %g with bandwidth %g and binsize %g ',sum(pdf_est).*binsize,bndw,binsize)
 k_earn2=ksdensity(earn,earnpoints);  
plot(earnpoints,[pdf_est, k_earn2]);


% Something important is happening around 0 and above 200K 
% What part of the observations are we missing?