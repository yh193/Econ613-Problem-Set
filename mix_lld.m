function neg_llh= mix_lld(beta_input)
% beta is the estimated input

% the first element in the first row is the beta for price
% all elements from the second colums to the end are beta for family income
% and others
beta1=beta_input(1,1);
beta2=[zeros(7,1),beta_input(:,2:end)];

% get the data from workspace
data1 = evalin('base', 'product10');
data2 = evalin('base', 'varbyproduct');
databyc_m=evalin('base', 'c_m');

% get the probability matrix
prob_matrix=data1*beta1+data2*beta2;

% therefore we write the loglikelihood 
loglh=sum(sum(log(exp(prob_matrix)./sum(exp(prob_matrix),2)).*databyc_m));
neg_llh=-loglh;
end