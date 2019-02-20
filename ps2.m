%% Exercise 1 Data creation
% a)X1 is the vector from a uniform distribution with range 1:3
X1=unifrnd(1,3,[10000,1]);

% b)X2: vector from a gamma distribution with shape 3 and scale 2
X2=gamrnd(3,2,[10000,1]);

% c)X3: vector from a binomial distribution with probability 0.3
X3=binornd(1,0.3,[10000,1]);

% d)eps: vector from a normal distribution with mean 2 and sd 1
eps=normrnd(2,1,[10000,1]);

% e)create the variable Y
Y=0.5+1.2*X1+(-0.9)*X2+0.1*X3+eps;

% f)create the variale ydum
mean_y=ones(10000,1);
mean_y=mean(Y)*mean_y;
ydum= Y>mean_y;

%% Exercise 2 OLS
% a)correlation between Y and X1
corr=corrcoef(Y,X1);
corr
% it is clearly to see the correlation betwen Y and X1 is far from 1.2

% b)the coefficients on the regression is 
X=[ones(10000,1),X1,X2,X3];
coeff_ols=(inv((X')*(X)))*(X')*(Y);
coeff_ols

% c)OLS standard error
% firstly calculating the sum of square error omega
omega=sum((Y-X*coeff_ols).^2);
% decondly using the formula for standard deviation of parameters
coeff_var=(inv((X')*X))*(X')*omega*X*(inv(X'*(X)));
coeff_std=sqrt(diag(coeff_var));

% d)standard error of bootstrap X
% datasample function set replacement as default
X_bootstrap=[X,Y];
  % d-1)bootstrap with 49 relications
  % get the vector to store the 49 estimators 
  coeff_bt1=ones(49,4);
  % first get the bootstrap sample with 10000 draws for 49 times
  for i=1:49
      data1=datasample(X_bootstrap,10000);
      x1=data1(:,1:4);
      y1=data1(:,5);
      coeff_olsbt1=(inv((x1')*(x1)))*(x1')*(y1);
      coeff_bt1(i,:)=coeff_olsbt1';
  end
  % get the std
  coeff_std1=std(coeff_bt1);
  coeff_std1

  % d-2)bootstrap with 499 relications
  % get the vector to store the 49 estimators
  coeff_bt2=ones(499,4);
  % first get the bootstrap sample with 10000 draws for 499 times
  for j=1:499
      data2=datasample(X_bootstrap,10000);
      x2=data2(:,1:4);
      y2=data2(:,5);
      coeff_olsbt2=(inv((x2')*(x2)))*(x2')*(y2);
      coeff_bt2(j,:)=coeff_olsbt2';
  end
  % get the std
  coeff_std2=std(coeff_bt2);
  coeff_std2

%% Exercise 3 Numerical Optimization
% a)see function that returns the likelihood of the probit in
% probit_likelihood.m

% b)steepest ascent optimization algorithm
% first, set the start point
beta0=[0.5;1.2;-0.9;0.1];
beta_ori=beta0;
% second, get the gradient using probit_likelihood.m
% set a vector to store the gradient
gradient=ones(4,1);

% if the difference is less than 0.000001, break the loop
diff=1;
while diff>0.000001 %steepest asent terminating condition
    
    for m=1:4
        % set h as a very tiny change
        h=0.01;
        beta1=beta_ori;
        beta2=beta_ori;
        beta1(m,1)=beta_ori(m,1)+h;
        beta2(m,1)=beta_ori(m,1)-h;
        gradient(m,1)=(probit_likelihood(X,ydum,beta1)-probit_likelihood(X,ydum,beta2))/(2*h);
    end
    %get the best scaling parameter in order to calculate best beta
    % set the starting point for optimizing best scaling parameter
    alpha=0;
    % write the likelihood funciton,
    ydumbar=ones(10000,1);
    ydumnot=(ydum<ydumbar);
    fun_best=@(x)-sum(ydum.* log(normcdf(X*(beta_ori-x.*gradient)))+ydumnot.*log(1-normcdf(X*(beta_ori-x.*gradient))));
    % change maximizing problem to minimizing problem by adding negative sign
    % get the optimization
    beta_probit_best=fminsearch(fun_best,alpha);
    % get the best beta
    beta_new=beta_ori-beta_probit_best.*gradient;
    % checking the difference of log likelihood between original and new beta
    diff=probit_likelihood(X,ydum,beta_new)-probit_likelihood(X,ydum,beta_ori);
    % assign the value for the next round steepest asent
    beta_ori=beta_new;
end

% showing the optimizaiton by steepest ascent
beta_new

% c)it is easily to see the distance between the estimator and ture parameter
distance=beta_new-beta0;
distance
% therefore, the parameter for intercept is far from the true parameter but
% the parameters for variables are close to the true parameters.

%% Exercise 4 Discrete choice
% get bootstrap datasample for ydum here
bt_ydum=[X,ydum];

% a)
% a-1)the probit model
% set the starting point
beta0=[0.5;1.2;-0.9;0.1];
% given the likelihood funciton we have in Exercise 3
ydumbar=ones(10000,1);
ydumnot=(ydum<ydumbar);
fun=@(x)-sum(ydum.* log(normcdf(X*x))+ydumnot.*log(1-normcdf(X*x)));
% change maximizing problem to minimizing problem by adding negative sign
% get the optimization  
beta_probit=fminsearch(fun,beta0);

% using 49 bootstrap to get the std of the parameter and calculate the t-stats
% set the vector to store the parameter value
bt_probit=ones(49,4);
% do bootstrap for 49 times
for i=1:49
%first get the bootstrap sample with 10000 draws 
data3=datasample(bt_ydum,10000);
x3=data3(:,1:4);
y3=data3(:,5);
% set the starting point
beta0=[0.5;1.2;-0.9;0.1];
% given the likelihood funciton we have in Exercise 3
ydumbar=ones(10000,1);
ydumnot1=(y3<ydumbar);
fun1=@(x)-sum(y3.* log(normcdf(x3*x))+ydumnot1.*log(1-normcdf(x3*x)));
% change maximizing problem to minimizing problem by adding negative sign
% get the optimization  
beta_probit1=fminsearch(fun1,beta0);
bt_probit(i,:)=beta_probit1';
end
% get the std 
probit_std1=std(bt_probit);
% get the tstas
probit_tstat=beta_probit./probit_std1';

% a-2)the logit model
% set the starting point
beta0=[0.5;1.2;-0.9;0.1];
% given the likelihood funciton we have in Exercise 3
ydumbar=ones(10000,1);
ydumnot=(ydum<ydumbar);
% creat logistic distribution object
logistic=makedist('logistic','mu',0,'sigma',1);
% write the log likelihood function using the logistic distribution object
fun_logit=@(x)-sum(ydum.* log(cdf(logistic,X*x))+ydumnot.*log(1-cdf(logistic,X*x)));
% change maximizing problem to minimizing problem by adding negative sign
% get the optimization  
beta_logit=fminsearch(fun_logit,beta0);

% using 49 bootstrap to get the std of the parameter and calculate the t-stats
% set the vector to store the parameter value
bt_logit=ones(49,4);
% do bootstrap for 49 times
for i=1:49
%first get the bootstrap sample with 10000 draws 
data4=datasample(bt_ydum,10000);
x4=data4(:,1:4);
y4=data4(:,5);
% set the starting point
beta0=[0.5;1.2;-0.9;0.1];
% given the likelihood funciton we have in Exercise 3
ydumbar=ones(10000,1);
ydumnot2=(y4<ydumbar);
% creat logistic distribution object
logistic=makedist('logistic','mu',0,'sigma',1);
% write the log likelihood function using the logistic distribution object
fun_logit=@(x)-sum(y4.* log(cdf(logistic,x4*x))+ydumnot2.*log(1-cdf(logistic,x4*x)));
% change maximizing problem to minimizing problem by adding negative sign
% get the optimization  
beta_logit1=fminsearch(fun_logit,beta0);
bt_logit(i,:)=beta_logit1';
end
% get the std 
logit_std=std(bt_logit);
% get the tstas
logit_tstat=beta_logit./logit_std';

% a-3)the linear probability model
% set the starting point
beta0=[0.5;1.2;-0.9;0.1];
% given the likelihood funciton we have in Exercise 3
ydumbar=ones(10000,1);
ydumnot=(ydum<ydumbar);
% write the log likelihood function
fun_linear=@(x)-sum(ydum.* log(X*x)+ydumnot.*log(1-X*x));
% change maximizing problem to minimizing problem by adding negative sign
% get the optimization  
beta_linear=fminsearch(fun_linear,beta0);

% using 49 bootstrap to get the std of the parameter and calculate the t-stats
% set the vector to store the parameter value
bt_linear=ones(49,4);
% do bootstrap for 49 times
for i=1:49
%first get the bootstrap sample with 10000 draws 
data5=datasample(bt_ydum,10000);
x5=data5(:,1:4);
y5=data5(:,5);
% set the starting point
beta0=[0.5;1.2;-0.9;0.1];
% given the likelihood funciton we have in Exercise 3
ydumbar=ones(10000,1);
ydumnot3=(y5<ydumbar);
% write the log likelihood function using the logistic distribution object
fun_linear=@(x)-sum(y5.* log(x5*x)+ydumnot3.*log(1-x5*x));
% change maximizing problem to minimizing problem by adding negative sign
% get the optimization  
beta_linear1=fminsearch(fun_linear,beta0);
bt_linear(i,:)=beta_linear1';
end
% get the std 
linear_std=std(bt_linear);
% get the tstas
linear_tstat=beta_linear./linear_std';

% b)compare the estimators (first colum) and their significances(second column)
probit_result=[beta_probit,probit_tstat];
logit_result=[beta_logit,logit_tstat];
linear_result=[beta_linear,linear_tstat];
probit_result
% Except the intercept, estimators of other three coefficients are close to
% real values and all of four estimators are significant, where the tstats
% are bigger than 2.
logit_result
% Estimators of four coefficients are not very close to real values but all 
% of four estimators are significant, where the tstats are bigger than 2.
linear_result
% Estimators of four coefficients are far from real values and all of 
% four estimators are not significant, where the tstats are far less than 2.


%% Exercise 5 Marginal Effects
% a)compute the marginal effect of X on Y for probit and logit models
% Given the formula of marginal effect, we have
  % a-1)for probit
  me_probit=normpdf(X*beta_probit)*beta_probit';
  % a-2)for logit
  me_logit=pdf(logistic,X*beta_logit)*beta_logit';
  
% b)computing the std of marginal effect
  % b-1) The delta method
  % Probit
  % get the variance-covariance matrix for beta
  % recall the 49 bootstrap beta in Exercise 4 
  Cov_probit=cov(bt_probit);
  % get the Jacobian matrix
  % first, get the Jacobian for each observation
    %get a very tiny change t
    t=0.000001;
    Jb_probit=ones(4,4);
    for p=1:4
        beta_me_probit1=beta_probit;
        beta_me_probit2=beta_probit;
        %get derivatives for each parameter
        beta_me_probit1(p)=beta_probit(p)+t;
        beta_me_probit2(p)=beta_probit(p)-t;
        %calculate each column of Jacobian
        Jb1=(normpdf(X*beta_me_probit1)*beta_me_probit1'-normpdf(X*beta_me_probit2)*beta_me_probit2')./(2*t);
        Jb_probit(:,p)=mean(Jb1)';
    end
  % given the formula of delta method for variance matrix of beta
  Var_probit=Jb_probit'*Cov_probit*Jb_probit;
  std_dm_probit=sqrt(diag(Var_probit));
  std_dm_probit
  
  % Logit
  % get the variance-covariance matrix for beta
  % recall the 49 bootstrap beta in Exercise 4 
  Cov_logit=cov(bt_logit);
  % get the Jacobian matrix
   % first, get the Jacobian for each observation
    %get a very tiny change t
    t=0.000001;
    Jb_logit=ones(4,4);
    for p=1:4
        beta_me_logit1=beta_logit;
        beta_me_logit2=beta_logit;
        %get derivatives for each parameter
        beta_me_logit1(p)=beta_logit(p)+t;
        beta_me_logit2(p)=beta_logit(p)-t;
        %calculate each column of Jacobian
        %pdf(logistic,X*beta_logit)*beta_logit'
        Jb2=(pdf(logistic,X*beta_me_logit1)*beta_me_logit1'-pdf(logistic,X*beta_me_logit2)*beta_me_logit2')./(2*t);
        Jb_logit(:,p)=mean(Jb2)';
    end
  % given the formula of delta method for variance matrix of beta
  Var_logit= Jb_logit'*Cov_logit* Jb_logit;
  std_dm_logit=sqrt(diag(Var_logit));
  std_dm_logit

  % b-2) Bootstrap
  %Probit
  %using 49 bootstrap to get average ME for 49 times
  % do bootstrap for 49 times
  for i=1:49
      %first get the bootstrap sample with 10000 draws
      data6=datasample(bt_ydum,10000);
      x6=data6(:,1:4);
      y6=data6(:,5);
      % set the starting point
      beta0=[0.5;1.2;-0.9;0.1];
      % given the likelihood funciton we have in Exercise 3
      ydumbar=ones(10000,1);
      ydumnot4=(y6<ydumbar);
      fun_probit_me=@(x)-sum(y6.* log(normcdf(x3*x))+ydumnot4.*log(1-normcdf(x3*x)));
      % change maximizing problem to minimizing problem by adding negative sign
      % get the optimization of estimators
      beta_probit_me=fminsearch(fun_probit_me,beta0);
      % get the ME
      me_probit1=normpdf(x6*beta_probit_me)*beta_probit_me';
      avg_me=mean(me_probit1);
      me_probit(i,:)=avg_me';
  end
  % get the std
  me_probit_std=std(me_probit);
  me_probit_std

  % Logit
  %using 49 bootstrap to get average ME for 49 times
  % do bootstrap for 49 times
  for i=1:49
      %first get the bootstrap sample with 10000 draws
      data7=datasample(bt_ydum,10000);
      x7=data7(:,1:4);
      y7=data7(:,5);
      % set the starting point
      beta0=[0.5;1.2;-0.9;0.1];
      % given the likelihood funciton we have in Exercise 3
      ydumbar=ones(10000,1);
      ydumnot5=(y7<ydumbar);
      fun_logit_me=@(x)-sum(y7.* log(cdf(logistic,x7*x))+ydumnot5.*log(1-cdf(logistic,x7*x)));
      % change maximizing problem to minimizing problem by adding negative sign
      % get the optimization of estimators
      beta_logit_me=fminsearch(fun_logit_me,beta0);
      % get the ME
      me_logit1=pdf(logistic,X*beta_logit_me)*beta_logit_me';
      avg_me_logit=mean(me_logit1);
      me_logit(i,:)=avg_me_logit';
  end
  % get the std
  me_logit_std=std(me_logit);
  me_logit_std



