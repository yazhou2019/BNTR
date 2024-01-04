# BNTR: Broadcasted Nonparametric Tensor Regression
- This package implements the method proposed in Y. Zhou, Wong, & He (2020).

  


# Installation
This package can be installed via function install_github in R package devtools:
```markdown
#install.packages("./BNTR_0.1.0.tar.gz", repos = NULL, type = "source")
install.packages("devtools")
devtools::install_github("yazhou2019/BNTR/BNTR")
```

# Experiments  
## Experiment 1: A toy example for a 32-by-32 input (not in the paper)
The following is a simple example for the 32-by-32 input. It can be implemented after the package BNTR is installed.

```markdown
library(BNTR)


#######################
#### simulate data ####
#######################


# load the true coefficient tensor (matrix)
data("X_horse")
BB <- X_horse
set.seed(2019)
# sample size of the training and test set
n_train <- 400
n_test <- 100
# signal level
signal_level = 0.5
# generate the training data
X_train = array(runif(prod(c(dim(BB), n_train)), 0, 1), c(dim(BB), n_train))
BroX_train = X_train + 0.6 * sin(2 * pi * (X_train - 0.5)^2)
y_train = 1 + crossprod(matrix(BroX_train, c(prod(dim(BB)), n_train)), as.vector(BB)) + signal_level * rnorm(n_train)
# generate the test data
X_test = array(runif(prod(c(dim(BB), n_test)), 0, 1), c(dim(BB), n_test))
BroX_test = X_test + 0.6 * sin(2 * pi * (X_test - 0.5)^2)
y_test = 1 + crossprod(matrix(BroX_test, c(prod(dim(BB)), n_test)), as.vector(BB)) + signal_level * rnorm(n_test)



###########################################
##### BroadcasTR ####
###########################################


############tuning paramters###############
# the tuning parameters have been selected by the hold-out method; see the paper
# rank used in the algorithm
rank <- 4
# the tuning parameters of the elastic net penalty
lambda_1 <- 5
lambda_2 <- 0.5  # 0 corresponds to lasso; 1 corresponds to ridge
# use the initial point from the sequential down-sizing procedure (1) or not (0). if not, we adopt a sequential warmstart descriped in the paper.
input_initial <- 0
if (input_initial == 1) {
 data("initial_point")
 beta0 <- initial_point$beta0
 B0 <- initial_point$B0
 warmstart <- 0
} else {
beta0 <- NA
 B0 <- NA
 warmstart <- 1
}
# number of basis
num_knots <- 5
knots = quantile(c(X_train), probs = c(seq(0, 1, 1/(num_knots - 1))))


###########training##################
tildePhiX_train = tildePhiX_trans(X_train, knots)
res = broadcasted_sparsetenreg(tildePhiX_train, y_train, r = rank, lambda = lambda_1, alpha = lambda_2, warmstart = warmstart, beta0 = beta0, B0=B0, Replicates=1)


###########important region##########
# norm tensor
normtensor <- fhatnorm_ten(full_R(res$beta), knots)
par(mar = c(0, 0, 0, 0), mfcol = c(1, 2), mai = c(0.2, 0.4, 0.4, 0.2))
plot(c(1, dim(BB)[1]), c(1, dim(BB)[2]), xaxt = "n", yaxt = "n", type = "n")
rasterImage(BB, 0, 0, dim(BB)[1] + 1, dim(BB)[2] + 1)
mtext("True coefficient tensor (matrix)", 3, line = 0.2)
plot(c(1, dim(BB)[1]), c(1, dim(BB)[2]), xaxt = "n", yaxt = "n", type = "n")
mtext("Estimated norm tensor (matrix)", 3, line = 0.2)
rasterImage(normtensor, 0, 0, dim(BB)[1] + 1, dim(BB)[2] + 1)


###########test#####################
# prediction on the test set
tildePhiX_test = tildePhiX_trans(X_test, knots)
y_pre = res$beta0 + crossprod(matrix(tildePhiX_test, c(prod(dim(full_R(res$beta))), n_test)), as.vector(full_R(res$beta)))
# prediction performance
cat("The prediction performance in these tunning parameters: \n", "MSPE =", sum((y_test - y_pre)^2)/n_test, "\n")


```

## Experiment 2: The synthetica data in Section 5.1 of the main paper
For the code to reproduce Table 1 in the main paper, please refer to ./experiments/Table_1_new. You can obtain the simulation results by using the following steps. 
#### For BroadcasTR:
- Step 1. Set ./experiments/Table_1_new/BrodcasTR as the working directory. 
- Step 2. Run the following code to generate the synthetic data. 
```markdown
nohup  Rscript 20230918_generate_data_new_editor.R. #You will obtain the synthetic data "Simul_n1000_rep50_final_fix_new_editor.Rdata" in "./SimResults".
```
- Step 3. Run the following code in the command line.
```markdown
nohup  Rscript --vanilla "SimNonLin.R" > ./logs 2>&1 & # obtain the fitting results
nohup  Rscript --vanilla "ISE_SimNonLin.R" > ./ISE_logs 2>&1 & # compute the ISE, after you obtain the fitting results
```
<details>
  <summary>Click to view collapsible paragraph---other methods </summary>
  
#### For TLR-rescaled 
- Step 1. Put the generated synthetic data in ./experiments/Table_1_new/TLR-rescaled/SimResults
- Step 2. Set ./experiments/Table_1_new/TLR-rescaled as the working directory.
- Step 3. Run the following code in the command line.
```markdown
nohup  Rscript --vanilla "SimNonLin.R" > ./logs 2>&1 &
nohup  Rscript --vanilla "ISE_SimNonLin.R" > ./ISE_logs 2>&1 &
```

#### For ENetR: 
- Step 1. Put the generated synthetic data in ./experiments/Table_1_new/ENetR /SimResults
- Step 2. Set ./experiments/Table_1_new/ENetR as the working directory.
- Step 3. Run the following code in the command line.
```markdown
nohup  Rscript --vanilla "SimENetR.R" > ./logs 2>&1 &
nohup  Rscript --vanilla "ISE_SimENetR.R" > ./ISE_logs 2>&1 &
```

#### For TLR: 
- Plase refer to the "TLR1" section in./experiments/Table_1_new/README.md.

 </details>
 
- Remark 1: You need to install the R package snow (for the parallel computing) and the dependences of BNTR in a linux server.
- Remark 2: You can set the tuning parameters of BroadcasTR in ./experiments/Table_1_new/ParallelComput/parallel_source1000.R. for the sample size n=1000.
- Remark 3: You can set the number of CPUs for the computation in ./experiments/Table_1_new/ParallelComput/parallel_replications_big_1000K8_new.R for the sample size n=1000. 
- Remark 4: For more details, please refer to ./experiments/Table_1_new.
- Remark 5: For Fig. 4 in the main paper, please refer to ./experiments/Figure_4_new. 
- Remark 6: For Fig. S.1. and Fig. S.2. in the supplementary material, please refer to ./experiments/Figure_5_new and ./experiments/Figure_6_new, respectively.


## Experiment 3: The real data analysis in Section 5.2 of the main paper
For the code to reproduce Table 2 in the main paper, please refer to ./experiments/Table 2. You can obtain the results by using the following steps. 
- Step 1. Download the data from http://neurotycho.org/data/20100802s1epiduralecogfoodtrackingbkentaroshimoda 
- Step 2. For data preprocessing, please follow the steps in ./experiments/Table 2/ECOG_observation/README.md. After you finish the data preprocessing, you will obtain X_test_1.Rdata, X_train_1.Rdata, and y_all.Rdata. Move them into the ./data directory.
- Step 3. Set ./experiments/Table 2 as the working directory. 
- Step 4. Run the following code in the command line.
```markdown
nohup  Rscript --vanilla "RealNonLin.R" > ./BroadcasTR_logs 2>&1 & # for BroadcasTR
nohup  Rscript --vanilla "RealTLR-rescaled.R" > ./RealTLR2_logs 2>&1 & # for TLR-rescaled
nohup  Rscript --vanilla "RealENet.R" > ./ENet_logs 2>&1 & # for ENetR
```
As for TLR, please refer to the "TLR1" section of ./experiments/Table 2/README.md. 

- Remark 1: You also need to install the R package snow (for the parallel computing) and the dependences of BNTR in a linux server.
- Remark 2: The data spliting follows ./data/idtest_matrix.Rdata and ./data/idtrain_matrix.Rdata. 
- Remark 3: For more details, please refer to ./experiments/Table 2.



## Experiment 4: The simulated monkey electrocorticography data in Section 5.3 of the main paper
For the code to reproduce Table 3 in the main paper, please refer to ./experiments/Table 3/RealDataSim. You can obtain the results by using the following steps. 
- Step 1. Move X_test_1.Rdata, X_train_1.Rdata, and y_all.Rdata (generated in the step 2 of Experiment 3) and monkey_1_new_K11_tuning.Rdata (generated in the step 4 of Experiment 3) into ./experiments/Table 3/RealDataSim/data. 
- Step 2. Set ./experiments/Table 3/RealDataSim as the working directory.
- Step 3. Run the following code to generate the simulated monkey electrocorticography data. 
```markdown
Rscript generate_real_sim_Data_new.R
```
- Step 4. Run the following code in the command line.
```markdown
nohup  Rscript --vanilla "RealNonLin_repro_K11_tuning.R" > ./BroadcasTR_logs 2>&1 & # for BroadcasTR
nohup  Rscript --vanilla "RealTLR-rescaled.R" > ./RealTLR2_logs 2>&1 & # for TLR-rescaled
nohup  Rscript --vanilla "RealENet.R" > ./ENet_logs 2>&1 & # for ENetR
```
As for TLR, please refer to the "TLR1" section of ./experiments/Table 3/RealDataSim/README.md. 

- Remark 1: You also need to install the R package snow (for the parallel computing) and the dependences of BNTR in a linux server.
- Remark 2: For more details, please refer to ./experiments/Table 3/RealDataSim.
- Remark 3: For Fig. S.3. in the supplementary material, please refer to ./experiments/Table 3. 

## Other Experiments : Extra simulations in the supplementary material
For other experiments, please refer to the "File Description List" secition. For example, the item (Table S.1. in the supplementary material: "./experiments/Table 1") implies that the code to reproduce Table S.1 is located in the "./experiments/Table 1" directory. 



# File Description List
In summary, all the code can be found in ./experiments, in which each subdirectory corresponds to figures or tables in the paper. 


### The main paper (the current version under review) 
- Table 1 in the main paper："./experiments/Table_1_new".
  
- Table 2 in the main paper: "./experiments/Table 2". 

- Table 3 in the main paper: "./experiments/Table 3".
  
- Fig. 4 in the main paper: "./experiments/Figure_4_new". 
  
### The supplementary material (the current version under review) 
<details>
  <summary>Click to view collapsible paragraph---The supplementary material (The latest)</summary>
 
- Fig. S.1. in the supplementary material: "./experiments/Figure_5_new".
  
- Fig. S.2. in the supplementary material: "./experiments/Figure_6_new". 

- Fig. S.3. in the supplementary material: "./experiments/Figure 7".

- Table S.1. in the supplementary material: "./experiments/Table 1".
  
- Table S.2 in the supplementary material: "./experiments/Supp/Table S1-S3".
  
- Fig. S.4. in the supplementary material: "./experiments/Figure 4".
  
- Fig. S.5. in the supplementary material: "./experiments/Figure 5".
  
- Fig. S.6. in the supplementary material: "./experiments/Figure 6".
  
- Fig. S.7. in the supplementary material: "./experiments/Supp/Figure S5". 

- Fig. S.8. in the supplementary material: "./experiments/Supp/Figure S6-S7". 

- Fig. S.9. in the supplementary material: "./experiments/Supp/Figure S6-S7".

</details>
 

# References
- Zhou, Y., Wong, R. K. W., & He, K. (2020). Broadcasted nonparametric tensor regression. arXiv preprint arXiv:2008.12927. [\[link\]](https://arxiv.org/abs/2008.12927v2)
- Remark: there is an updated version which is under review.
