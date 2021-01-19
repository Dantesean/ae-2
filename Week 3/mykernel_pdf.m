function [k_pdf_hat] = mykernel_pdf(X,X0,h)

% INPUTS
    %  X          : Nx1  vector of data  (if X contains more columns only the
    %  X0         : point where density is to be evaluated
    %  h          : bandwidth
%  (c) 
  
% OUTPUTS
    % k_pdf_hat   : kernel Density estimate at X0
% local variables used
    % Z , nobs
    [nobs,ncols]  = size(X);
    Z = (X - X0)/h;   % look up formula
    kZ = max((3/4*(1-Z.^2)), zeros(nobs,1));
    k_pdf_hat = mean(kZ)/h;   % look up formula and calculate the pdf^  Hint: use the mean of K(Z) and divide h
end