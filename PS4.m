%%Problem Set 4
%%Yaqi Hu

% deal with the format of the data
data=table2array(KoopTobias);

%% Exercise 1
% Selecte the ID from 1 to 5 
id=data(:,1)<6;
data1to5=data(id,:);
% Take the variable EDUC(Education) which is in the second column as an example
data_edu=data1to5(:,1:2);
data_edu
% Therefore it is easily to see that this is an unbalanced data. 

%% Exercise 2
% get the independent variables 
X1=[data(:,2),data(:,4)];
Y1=data(:,3);

% First stage regression
% run the OLS by function fitlim.m
mdl1=fitlm(X1,Y1);
% get the variance of the estimated ramdom effect alpha
std_alpha=table2array(mdl1.Coefficients(1,2));
var_alpha=std_alpha^2;

% get the 17919 raw residuals
res1=table2array(mdl1.Residuals(:,1));
% calculate the variance of it 
var_res1=var(res1);

% calculate the lamda to conduct the second step regression
% first, calculate the maximum number of time period
T=data(:,1);
T_count = hist(T,unique(T));
T_max=max(T_count);
lamda=1-(var_res1/(var_res1+T_max*var_alpha));

% Second step regression
% Deal with the data first
% get original data 
data1=[data(:,1),data(:,2:4)];
% set the matrix to contain the data after calculation
data2=[data(:,1),zeros(17919,3)];

for i=1:2178
    id1=find(data1(:,1)==i);
    data_inter=data1(id1,:);
    data_mean=mean(data_inter);
    data_final=data_inter-lamda*data_mean;
    data2(id1,2:4)=data_final(:,2:4);  
end

% get the data we need
X2=[data2(:,2),data2(:,4)];
Y2=data2(:,3);
% estimate by OLS
mdl2=fitlm(X2,Y2);
% get the estimated beta and t statistics
mdl2.Coefficients(2:3,1:3)

%%% for GLS estimation, we can use fucntion as well
[coeff,se] = fgls(X1,Y1);
[coeff,se]

%% Exercise 3
% Between Estimator
% Deal with the data first
% set the matrix to contain the data after calculation
data3=zeros(2178,3);
for i=1:2178
    id1=find(data1(:,1)==i);
    data_inter=data1(id1,:);
    data_mean=mean(data_inter,1);
    data3(i,:)=data_mean(1,2:4);  
end

% get the data we need
X3=[data3(:,1),data3(:,3)];
Y3=data3(:,2);
% estimate by OLS
mdl3=fitlm(X3,Y3);
% get the estimated beta and t statistics
mdl3.Coefficients(2:3,1:3)

% Within Estimator
% Deal with the data first
% set the matrix to contain the data after calculation
data4=[data(:,1),zeros(17919,3)];
for i=1:2178
    id1=find(data1(:,1)==i);
    data_inter=data1(id1,:);
    data_mean=mean(data_inter);
    data_final=data_inter-data_mean;
    data4(id1,2:4)=data_final(:,2:4);  
end

% get the data we need
X4=[data4(:,2),data4(:,4)];
Y4=data4(:,3);
% estimate by OLS
mdl4=fitlm(X4,Y4);
% get the estimated beta and t statistics
mdl4.Coefficients(2:3,1:3)

% First difference Estimater
% Deal with the data first
cell_data5=cell(2178,1);
for i=1:2178
data_inter=data1(find(data1(:,1)==i),:);
size_data=size(data_inter);

% set the matrix to store the difference data
if size_data(1)>1
    rows=size_data(1)-1;
data_diff=ones(rows,size_data(2));
 for j=1:rows
     data_diff(j,:)=data1(j+1,:)-data1(j,:);
 end

else
data_diff=data_inter;
end
% the cell to contain the data after calculation
cell_data5{i}=data_diff;

end

% change the data from the cell to matrix
data5=cell2mat(cell_data5);

% get the data we need
X5=[data5(:,1),data5(:,3)];
Y5=data5(:,2);
% estimate by OLS
mdl5=fitlm(X5,Y5);
% get the estimated beta and t statistics
mdl5.Coefficients(2:3,1:3)

%% Exercise 4
% set matrix to reflect individual fixed effect,it should have 2178 columns
% for original data
alpha=zeros(17919,2178);
for i=1:2178
    id1=find(data(:,1)==i);
    alpha(id1,i)=1; 
end
% combine the parameter with original data and keep 2 variant variables and
% all invariant variables
data6=[data(:,1),alpha,data(:,2),data(:,4),data(:,6:10),data(:,3)];

%% Part one
% select 100 ID data from ID 1 to ID 100
id4=data(:,1)<101;
data7=data6(id4,:);

% First step- OLS 
% get the X and Y
X7=data7(:,2:2181);
Y7=data7(:,2187);
mdl7=fitlm(X7,Y7);

% get the estimated fixed effet parmeter
alpha_est=table2array(mdl7.Coefficients(2:101,1));

% Second step- run regression on fixed effect by invariant variables
% get the invariant variables
var_inv=data7(:,2182:2186);
X7_inv=unique(var_inv,'stable','rows');
% regress it by function
mdl7_alpha=fitlm(X7_inv,alpha_est);
mdl7_alpha.Coefficients(2:6,1:3)

% Why it is not correctly estimated?
%the reason is that multistage estimation for individual fixed effect alpha
%contains the errors from the first step and therefore the output of std is
%not corret under formula.

%% Part two
% bootstrap the previos procedure for 20 times
% set the parameter to store the value 
beta_bootstrap=zeros(20,5);

for b=1:20
% select 100 ID data randomly from 2178
id_100 = randperm(2178,100);
data8=data6(ismember(data1(:,1),id_100),:);

% First step- OLS 
% get the X and Y
X8=data8(:,2:2181);
Y8=data8(:,2187);
mdl8=fitlm(X8,Y8);

% get the estimated fixed effet parmeter
alpha_all=table2array(mdl8.Coefficients(2:2179,1));
alpha_bootstrap=alpha_all(id_100,1);

% Second step- run regression on fixed effect by invariant variables
% get the invariant variables
var_bootstrap=[data8(:,1),data8(:,2182:2186)];
X8_inv=unique(var_bootstrap,'stable','rows');
% regress it by function
mdl8_alpha=fitlm(X8_inv,alpha_bootstrap);
beta_invariant=table2array(mdl8_alpha.Coefficients(2:6,1));

beta_bootstrap(b,:)=beta_invariant';
end

% catculating the std of five beta
std_beta=std(beta_bootstrap);
std_beta
