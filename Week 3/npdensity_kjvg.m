
function [Xmidpoints_used,k_pdf,bandwidth_used] = npdensity_kjvg(X,Xmidpoints,bandwidth)

% INPUTS
    %  X          : Nx1  vector of data  (if X contains more columns only the
    %                                      first column will be used)
    %  Xmidpoints : nrbins x 1 columnvector of bin midpoints
    %  bandwidth  :  if > 0 then this is the bandwidth that is used in the kernel density estimation   
    %                if 0 or <0 , let the program determine the bandwidth based
    %                on a plug-in estimator
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
    end;
    
    firstb = Xmidpoints_used(1,:);
    lastb  = Xmidpoints_used(nrbins,:);
    binsize = ((lastb - firstb)/nrbins); % average binsize in this case

if bandwidth>0; % use the bandwidth given in the proc argument
  bandwidth_used = bandwidth;
else 
  % SPECIFY THE BANDWITDTH(S)
  delta =0.7764;               % see (9.11) and table 9.1 Cameron & Trivedi
  bandwidth_used = 1.3643*delta*stdX*nr^(-1/5);
  % *(lastb-firstb)
  % use normal
  % use Silverman
  % use Cross Validation
end;    % endif 

k_pdf  = zeros(nrbins,1);  %  this will contain the density estimates at the midpoints_used
                           %  bin by bin calculations since reduces the amount of workspace used vs k_pdf at each observation %
  
    for J=1:nrbins;        % for each bin       %
       Xb = Xmidpoints_used(J,1);          % one bin at a time  %
       Z = (iota*Xb - X)/bandwidth_used;
       KX = pdf('Normal',Z,0,1)/bandwidth_used;     % CHANGE if YOU WANT A DIFFERENT KERNEL  %
       k_pdf(J,1) = mean(KX);
    end    % for; 

end
