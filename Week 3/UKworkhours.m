% read the data
% column:
% 1: hours worked
% 2: Age
% 3: Gender (1 if Female)
clear;
UKworkhrs =importdata('UKworkHrsAgeGender.txt');
yy = UKworkhrs.data;
hrs =yy(:,1);    % hours worked
hrs1mid=linspace(0,96,100)';

[hrspoints,pdf_est,bndw] = npdensity_kjvg(hrs,hrs1mid,1.05);

plot(hrspoints,pdf_est);
% check that the density integrates to 1
binsize = (hrspoints(2)-hrspoints(1)); % since equally spaced 
sprintf('The density integrates to %g with bandwidth %g and binsize %g ',sum(pdf_est).*binsize,bndw,binsize)

% Compare with Matlab's own density estimator
k_hrs2=ksdensity(hrs,hrspoints);  
plot(hrspoints,[pdf_est, k_hrs2]);


% do a nonparametric regression of hrs on age
age= yy(:,2);
age1mid=linspace(15,69,100)';  % generate 100 midpoints where the regression is to be evaluated
  [Xpoints,m_reg,bndw] = npregress_kjvg(hrs,age,age1mid,1.05); % !? why does this not make sense?!
 plot(Xpoints,m_reg);

