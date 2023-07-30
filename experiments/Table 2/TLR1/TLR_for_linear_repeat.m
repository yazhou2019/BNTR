%Hua Zhou's code 
%tensor  linear regression

clear;
% reset random seed
s = RandStream('mt19937ar','Seed',2);
RandStream.setGlobalStream(s);


format long;


X_train=load('X_train_1_MATLAB.mat');
X_train=X_train.X_train;
X_test=load('X_test_1_MATLAB.mat');
X_test=X_test.X_test;


X_all=zeros(64,10,10,10000);
X_all(:,:,:,1:4000)=X_train;
X_all(:,:,:,4001:10000)=X_test;


y_train=load('y_train_1_MATLAB.mat');
y_train=y_train.y_train;
y_test=load('y_test_1_MATLAB.mat');
y_test=y_test.y_test;


y_all=vertcat(y_train,y_test);

idtrain_matrix=load('idtrain_matrix_ECOG_MATLAB');
idtrain_matrix=idtrain_matrix.idtrain_matrix;
idtest_matrix=load('idtest_matrix_ECOG_MATLAB');
idtest_matrix=idtest_matrix.idtest_matrix;

%iter=1;
MSE_grid_val_pre=zeros(21,4,2,20);

for iter=1:20


X_train=X_all(:,:,:,idtrain_matrix(iter,:));
X_test=X_all(:,:,:,idtest_matrix(iter,:));
y_train=y_all(idtrain_matrix(iter,:));
y_test=y_all(idtest_matrix(iter,:));



lambda = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1,2.5,5,7.5,10,25,50,75,100,250,500,750,1000];
R=[1,2,3,4,5,6,7,8];
pentype = 'enet';
penparam = 1;

[~,L]=size(lambda);
[~,RR]=size(R);

MSE_validation=zeros(L*RR,1);
MSE_pre=zeros(L*RR,1);
MSE_grid_pre=zeros(L,RR);
MSE_grid_val=zeros(L,RR);


[~,~,~,n]=size(X_train);
X=ones(n,1);

for r=1:RR
    for l=1:L
%[~,beta_rk3,glmstats3,dev3] = kruskal_reg(X,tensor(X_train),y_train,R(r),'normal');

tic;
[beta0,beta_rk1,~,~] = kruskal_sparsereg(X,tensor(X_train),y_train,R(r),'normal',...
    lambda(l),pentype,penparam);
toc;

BB=double(beta_rk1);

y=zeros(6000,1);
for i=1:6000
y(i) = beta0 +(ttt(tensor(X_test(:,:,:,i)), tensor(BB), 1:3));
end

MSE_validation(L*(r-1)+l)=sum((y(1:1000)-y_test(1:1000)).^2)/1000;
MSE_pre(L*(r-1)+l)=sum((y(1001:6000)-y_test(1001:6000)).^2)/5000;

disp([R(r),lambda(l),MSE_validation(L*(r-1)+l),MSE_pre(L*(r-1)+l)]);

MSE_grid_val(l,r)=MSE_validation(L*(r-1)+l);
MSE_grid_pre(l,r)=MSE_pre(L*(r-1)+l);

    end
end

MSE_grid_val_pre(:,:,1,iter)=MSE_grid_val;
MSE_grid_val_pre(:,:,2,iter)=MSE_grid_pre;


end
%save('MSE_grid_val_pre_l1','MSE_grid_val_pre');
%save('MSE_grid_val_pre_l2','MSE_grid_val_pre');

%save('MSE_grid_val_pre_linear_10_l1','MSE_grid_val_pre');



[x,y]= find(MSE_grid_val_pre(:,:,1)==min(min(MSE_grid_val_pre(:,:,1))));
lambda(x)
R(y)
MSE_grid_val_pre(x,y,2)


MSE_all=zeros(20,1);
for iter=1:20
[x,y]= find(MSE_grid_val_pre(:,:,1,iter)==min(min(MSE_grid_val_pre(:,:,1,iter))));    
  disp([x(1),y(1)]) 
  MSE_all(iter,1)= MSE_grid_val_pre(x(1),y(1),2,iter);
  %disp(MSE_grid_val_pre(x(1),y(1),2,iter))
end

mean(MSE_all)
sqrt(var(MSE_all))
    
%save('MSE_linear_l1_20','MSE_all');

