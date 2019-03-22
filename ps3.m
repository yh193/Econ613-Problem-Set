%% Problem Set 3
% Name: Yaqi Hu 
% Net ID: yh193

% download data from the github and import the data
% change the table type data to numbers
demos=table2array(demos);
product=table2array(product);

%% Exercise 1
% Part 1 Average and dispersion in product characteristics
% which means that the mean and standard deviation of 10 products 

% get the price of 10 product
product10=product(:,4:13);

% get the mean 
product_mean=mean(product10);
% get the standard deviaiton
product_std=std(product10);

% Part 2 Market share and share by product characteristics
% which means calculating the frequancy of each product

% get the raw data of choices
product_choice=product(:,3);
% calculting the frequency
 % set a variable for storing the frequency
 product_share=zeros(10,1);
 % get the initial data
 product_choice_bar=product_choice;
 % get the frequency of ten products one by one
 for i=1:10
     product_fre= (i-1)*ones(4470,1)<product_choice_bar& product_choice_bar<(i+1)*ones(4470,1);
     product_share(i,1)=sum(product_fre);
 end 
 
% Part 3 Mapping between observed attributes and choices
% which means we need to find the relationship between the choices of the products and
% income

% first merge the data by the same income level
 % get the all income level
 income=demos(:,3);
 uni_income=unique(income); % we get 14 income level
 % get the id data by income level
 id_income=demos(:,2:3);
 % get the choice with the ID 
 product_byincome=product(:,2:3);
 
 % set the variable to store the choices for each income level
 choicebyincome=zeros(14,10);
 
  for j=1:14
      % get the all id respect to the specifc income level
      idex=find(id_income(:,2)==uni_income(j));
      id=id_income(idex);
    
      % get the all choices respect to the specifc id
      id_size=size(id);
      k=id_size(1);
      % store the product choices by idd for each income level
      choice_byincome=[];
      for l=1:k
      idex1=find(product_byincome(:,1)==id(k));
      productbyincome=product_byincome(idex1,2);
      choice_byincome=[choice_byincome;productbyincome];
      end
      
      % calculting the frequency as in part 2
         % set a variable for storing the frequency
         income_share=zeros(10,1);
         % get the initial data
         product_choice_bar=choice_byincome;
         size_choice=size(product_choice_bar);
         num=size_choice(1);
         % get the frequency of ten products one by one
         for i=1:10
           product_fre= (i-1)*ones(num,1)<product_choice_bar& product_choice_bar<(i+1)*ones(num,1);
           income_share(i,1)=sum(product_fre);
         end
         
       %store the market share of products by specific income level 
       choicebyincome(j,:)=income_share';        
  end

%% Exercise 2
% first model is the conditional logit model 

% writing the choice matrix
c_m=zeros(4470,10);
for n=1:4470
c_m(n,product_choice(n,1))=1;
end

% from the result of running the optimization function, we need to increase
% the MaxFunEvals and MaxIter option with default value 200*numberOfVariables. 
options = optimset('MaxFunEvals',200000,'MaxIter',200000);

% writing the log likelihood function
% deal with the data
data=product10-product10(:,1);
% assume the price effect parameter and product effect parameter 
fun_clogit=@(x)-sum(sum(log(exp(data*x(1)+[0,x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)])./sum(exp((data*x(1)+[0,x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)])),2)).*c_m));

% optimizing the model
% set the starting point
start=[0,0,0,0,0,0,0,0,0,0];
% set estimated variables
para_best=fminsearch(fun_clogit,start,options);

% Comment on the coefficient of the price
% Only the sign of the estimated parameter beta matters here and the magnitude
% is not important here.
% If the sign of beta is positive, it means the rise of the price will
% improve the probability that people choose any kind of product.


%% Exercise 3
% the second model is the multinomial model

% write the multi_lld.m function

% optimizing the model
options = optimset('MaxFunEvals',200000,'MaxIter',200000);

% set the starting point
parabystart=zeros(7,9);
% set estimated variables
para_ps3=fminsearch(@multi_lld,parabystart,options);

% Comment on the coefficient on family
% Only the sign of the estimated parameter family income matters here and the magnitude
% is not important here.
% If the sign of paramters is positive, it means the rise of the family income will
% improve the probability that people choose certain product as each
% product has an respective income parameter.

%% Exercise 4

% Part 1 marginal effect of conditional logit model
% given the formula from the slide, we have
 beta=para_best(1,1);
 alpha=[0,para_best(1,2:end)];
% get the probability maatrix first
 prob_esti=exp(data*beta+alpha)./sum(exp(data*beta+alpha),2);
 
 me_condi=(ones(10,1)*(sum(prob_esti)*beta)).*diag(ones(1,10))-(prob_esti'*prob_esti)*beta;
 mebycondi=me_condi/4470;
% Part 2 marginal effect of multinominal logit model
% given the formula fromt the slide, we have
 % get the original data
 var1=[demos(:,2:5),demos(:,7:9),ones(516,1)];
 var2=product(:,2:3);

 varbyproduct=zeros(4470,7);

for p=1:4470
    idexbydemos=find(var1(:,1)==var2(p,1));
    varbyproduct(p,:)=var1(idexbydemos,2:8);
end
 % firstly get the probability matrix
 beta_me=[zeros(7,1),para_ps3];
 p_esti=exp(varbyproduct*beta_me)./sum(exp(varbyproduct*beta_me),2);
 
 inter_me=p_esti*beta_me';
 me_multi=p_esti.*(ones(4470,1)*beta_me(1,:)-inter_me(:,1)*ones(1,10));
 % get the marginal effect of the income variables to different products
 mebymulti=mean(me_multi);

 % Comment on marginal effect
 % Not only the sign of marginal effect matters here but also the
 % maginitude of parameter matters. 
 % For the first model, it shows the change of probability people choose
 % each product if the price of any product changes;
 % For the second model, it shows the change of probability people choose
 % different products if the family income changes;
 
 
 %% Exercise 5
 % Part 1 optimize the mixed logit
 % write the mix_lld.m function

 % optimizing the model
 options = optimset('MaxFunEvals',200000,'MaxIter',200000);

 % set the starting point
 start_beta=[1,0,0,0,0,0,0]';
 start1=[start_beta,zeros(7,9)];
 % set estimated variables
 para_ps5_1=fminsearch(@mix_lld,start1,options);
 
 % Part 2 optimize the restricted mixed logit
 % write the mixre_lld.m function
 
 % optimizing the model
 options = optimset('MaxFunEvals',200000,'MaxIter',200000);

 % set the starting point
 start_beta=[1,0,0,0,0,0,0]';
 start2=[start_beta,zeros(7,8)];% deleting the sample who choose product 1
 % set estimated variables
 para_ps5_2=fminsearch(@mixre_lld,start2,options);
 
 % Part 3 Compute the statistics
 MTT=(-2)*(-mix_lld(para_ps5_1))-(-mixre_lld(para_ps5_2));
 
 % kai distribution is Gamma(U/2,1/2), the freedom here is 56
 prob_critical=gamcdf(MTT,28,0.5);
 
 % Comment on IIA
 % From the result of the MTT statistics, we reject the hypothesis that the
 % restricted model has the same beta as original beta and in other words,
 % IIA does not hold.
 
 
 
