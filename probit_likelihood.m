function [lh]=probit_likelihood(x,y,beta)
% x as the independant variables should have intercept
% y is the dependant variables and beta is the coefficient
% the number of columns of the x is the same as the numeber of 
% rows of the beta

% output lh is the log likelihood

% given the log likelihood function of probit model
% get the probability
p=normcdf(x*beta);
% get (1-y) term as ynot
s=size(y);
nrow=s(1);
ybar=ones(nrow,1);
ynot=(y<ybar);
% get each log likelihood
lh1=y.* log(p)+ynot.*log(1-p);
% get overall log likelihood
lh=sum(lh1);

end