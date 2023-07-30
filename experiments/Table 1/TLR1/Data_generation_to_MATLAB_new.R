library(R.matlab)


load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/Simul_n1000_rep50_final_fix_new.Rdata")

n = 1000
X_data_all <- array(0, c(64,64,n,50))
y_all_all <- array(0, c(4,n,50))

for(iter in 1:50){
X_data_all[,,,iter]=data_all[[iter]]$X
for (CASE in 1:5){
y_all_all[,,iter]=data_all[[iter]]$y[[CASE]] 
}
}





writeMat("X_data_all_mat.mat", X_data_all=X_data_all)
writeMat("y_all_all_mat.mat", y_all_all=y_all_all)



