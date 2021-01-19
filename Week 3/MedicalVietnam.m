% Medical expenditure in Vietnam
% Import the data on medical expenses and income in Vietnam
MedViet =importdata('Vietnxy.asc');
xy=MedViet.data;
% estimate the marginal densities
% run a non-parametric regression
