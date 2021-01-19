clear;
stockmr = importdata('STOCKMRCAPM.asc');
yx = stockmr.data;
retcompany = yx(:,1); % stock return of company A
retmarket = yx(:,2); % market stock return
[pdf] = mykernel_pdf(retcompany, retmarket, 2)