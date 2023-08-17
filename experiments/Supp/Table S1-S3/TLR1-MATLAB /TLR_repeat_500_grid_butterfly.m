%simulation of TLR
%clc;
%parpool('local', 2);








%rank
R=[1,2,3,4,5];
%tuning parameter 1
lambda=[0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000];
%tuning parameter 2
alpha=[0,0.5,1]+1;

MSE_all=zeros(100,5);
BB_all=zeros(64,64,5,50);
b0_all=zeros(5,50);


%clear X_data_all;
%clear X_data_train;
%clear X_data_test;
%clear X_data_test_matrix;
%clear X_data_train_matrix;

%parfor iter=1:50
    %pathdef;
for iter=1:50
    
X_data_all=load('X_data_all_mat.mat');
X_data_all=X_data_all.X_data_all;

y_data_all=load('y_all_all_mat.mat');
y_data_all=y_data_all.y_all_all;

X_data=X_data_all(:,:,:,iter);
y_data=y_data_all(:,:,iter);

clear X_data_all;
clear y_data_all;

n_use=500;
n_train=n_use * 0.8;
n_vali=n_use * 0.2;

y_train=y_data(:,1:n_train);
y_vali=y_data(:,((1000-n_vali+1):1000));
y_test=y_vali;

res_middle=cell(6);
res=cell(4);


%for signal_i =5
signal_i = 5
[res_middle{1},res_middle{2},res_middle{3},res_middle{4},res_middle{5},res_middle{6}]=validation_result_MATLAB_grid(R,alpha,lambda,X_data(:,:,1:n_train),y_train(signal_i,:),X_data(:,:,((1000-n_vali+1):1000)),y_vali(signal_i,:),X_data(:,:,((1000-n_vali+1):1000)),y_test(signal_i,:));
    res{signal_i}=res_middle;
    MSE_all(iter,signal_i)=res_middle{3}(1);
%BB_all{(signal_i-1)*50+iter}=res_middle{1};
%b0_all{(signal_i-1)*50+iter}=res_middle{6};
BB_all(:,:,signal_i,iter)=res_middle{1};
b0_all(signal_i,iter)=res_middle{6};

%end

end

csvwrite('MSE_all_500_grid_butterfly.csv',MSE_all)
save('BB_all_500_grid_butterfly.mat','BB_all')
save('b0_all_500_grid_butterfly.mat','b0_all')

%delete(gcp('nocreate'))
