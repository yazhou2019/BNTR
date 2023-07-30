%validation_result_MATLAB

function [BB,MSE_vali,MSE_pre,BBbig_nonlinear,b_validation_test_lambda_R_nonlinear,b0]=validation_result_MATLAB_grid(R,alpha,lambda,X_data_train,y_train,X_data_vali,y_vali,X_data_test,y_test)

%R=[1,3,5,7];
pentype = 'enet';
%penparam = 1;

num_R=length(R);
num_lambda=length(lambda);
num_alpha=length(alpha);
num_all=num_R*num_lambda*num_alpha;

[~,~,~, n]=size(X_data_train);
X=ones(n,1);


b_validation_test_lambda_R_nonlinear=zeros(num_all,6);
BBbig_nonlinear=cell(num_all);
for i1 = 1:num_R
for i2 = 1:num_alpha
for i3 = 1:num_lambda
if i3==1
[b0,beta,~,~] = kruskal_sparsereg(X,X_data_train,y_train,R(i1),'normal',...
    lambda(i3),pentype,alpha(i2));
else
[b0,beta,~,~] = kruskal_sparsereg(X,X_data_train,y_train,R(i1),'normal',...
                                  lambda(i3),pentype,alpha(i2),'B0',beta);
end


num_location=(i1-1)*num_alpha*num_lambda+(i2-1)*num_lambda+i3;
BB=double(beta); 

n_vali=length(y_vali);

y=zeros(n_vali,1);
for i=1:n_vali
y(i) = b0 +(ttt(tensor(X_data_vali(:,:,:,i)), tensor(BB), 1:3));
end

%size(y_vali)
%size(y)
MSE_vali=sum((y_vali-y).^2)/n_vali;
%size(MSE_vali)
num_location/num_all
n_test=length(y_test);

y=zeros(n_test,1);
for i=1:n_test
              y(i) = b0 +(ttt(tensor(X_data_test(:,:,:,i)), tensor(BB), 1:3));
end

MSE_test=sum((y_test-y).^2)/n_test;

BBbig_nonlinear{num_location}=BB;
b_validation_test_lambda_R_nonlinear(num_location,1)=b0;
b_validation_test_lambda_R_nonlinear(num_location,2)=MSE_vali;
b_validation_test_lambda_R_nonlinear(num_location,3)=MSE_test;
b_validation_test_lambda_R_nonlinear(num_location,4)=R(i1); 
b_validation_test_lambda_R_nonlinear(num_location,5)=alpha(i2);
b_validation_test_lambda_R_nonlinear(num_location,6)=lambda(i3);

    
end
end
end

[index_best]= find(b_validation_test_lambda_R_nonlinear(:,2)==min(b_validation_test_lambda_R_nonlinear(:,2)));
BB=BBbig_nonlinear{index_best(1)};
MSE_vali=b_validation_test_lambda_R_nonlinear(index_best(1),2);
MSE_pre=b_validation_test_lambda_R_nonlinear(index_best(1),3);
b0=b_validation_test_lambda_R_nonlinear(index_best(1),1);
end

