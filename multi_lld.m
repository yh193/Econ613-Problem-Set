function neg_llh=multi_lld(beta_input)
% the input is the beta parameters we want to estimate

% get the variable data from the workspace
databydemos = evalin('base', 'demos');
databyproduct= evalin('base', 'product');
databyc_m=evalin('base', 'c_m');

% get the six variables related to family income and adding the product
% effect variable
% var is a 516x7 matrix
var1=[databydemos(:,2:5),databydemos(:,7:9),ones(516,1)];
var2=databyproduct(:,2:3);
% we connect the information in demos with the product to expand the data
% sample number into 4470
% set an variable to store data
varbyproduct=zeros(4470,7);

for p=1:4470
    idexbydemos=find(var1(:,1)==var2(p,1));
    varbyproduct(p,:)=var1(idexbydemos,2:8);
end
% the beta_input is the 7x9 matrix we want to estimate
% beta is a 7x10 matrix with the first column all zeros
beta=[zeros(7,1),beta_input];

% therefore we write the loglikelihood 
loglh=sum(sum(log(exp(varbyproduct*beta)./sum(exp(varbyproduct*beta),2)).*databyc_m));
neg_llh=-loglh;

end
