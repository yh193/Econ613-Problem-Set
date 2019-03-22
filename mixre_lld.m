function neg_llh= mixre_lld(beta_input)
% beta is the estimated input

% the first element in the first row is the beta for price
% all elements from the second colums to the end are beta for family income
% and others
beta1=beta_input(1,1);
beta2=[zeros(7,1),beta_input(:,2:end)];

% get the data from workspace
choice= evalin('base', 'product_choice');
data1 = evalin('base', 'product10(:,2:end)');
data2 = evalin('base', 'varbyproduct');
databyc_m=evalin('base', 'c_m(:,2:end)');

% deletinng the samples who choose 1
datapool=[choice,data1,data2,databyc_m];

idex_product1 = find(datapool(:,1) == 1);  
datapool(idex_product1,:)=[];

data3=datapool(:,2:10);
data4=datapool(:,11:17);
data_c_m=datapool(:,18:26);
% get the probability matrix
prob_matrix=data3*beta1+data4*beta2;

% therefore we write the loglikelihood 
loglh=sum(sum(log(exp(prob_matrix)./sum(exp(prob_matrix),2)).*data_c_m));
neg_llh=-loglh;
end