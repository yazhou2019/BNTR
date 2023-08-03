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



y_all=load('y_all_sim2.mat');
y_all=y_all.y_all;

load('idtrain_matrix_MATLAB.mat')
load('idtest_matrix_MATLAB.mat')


MSE_all=zeros(10,2);
BB_all=zeros(64,10,10,10);
b0_all=zeros(1,10);

%rank
R=[1,2,3,4,5,6,7,8];
%tuning parameter 1
lambda=[0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1,2.5,5,7.5,10,25,50,75,100,250,500,750,1000];
%tuning parameter 2
alpha=[0,0.5,1]+1;




for iter=1:10
    pathdef;
 

res_middle=cell(6);
res=cell(4);


id_train = idtrain_matrix(iter, 1:4000);
id_vali = idtest_matrix(iter, 1:1000);
id_test = idtest_matrix(iter, 1001:6000);
%[BB,MSE_vali,MSE_pre,BBbig_nonlinear,b_validation_test_lambda_R_nonlinear,b0]

[res_middle{1},res_middle{2},res_middle{3},res_middle{4},res_middle{5},res_middle{6}]=validation_result_MATLAB_grid(R,alpha,lambda,X_all(:,:,:,id_train),y_all(id_train),X_all(:,:,:,id_vali),y_all(id_vali),X_all(:,:,:,id_test),y_all(id_test));
  

MSE_all(iter,1)=res_middle{2}(1);
MSE_all(iter,2)=res_middle{3}(1);
BB_all(:,:,:,iter)=res_middle{1};
b0_all(1,iter)=res_middle{6};


end

csvwrite('MSE_all_RealDataSim2.csv',MSE_all)
save('BB_all_RealDataSim2.mat','BB_all')
save('b0_all_RealDataSim2.mat','b0_all')

%delete(gcp('nocreate'))



