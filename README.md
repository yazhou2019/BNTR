# BNTR: Broadcasted Nonparametric Tensor Regression
- This package implements the method proposed in Y. Zhou, Wong, & He (2020).
- If you want to see the code for the analyses presented in the paper, please refer to the 'Experiments in the paper' section (In summary, all the code can be found in ./experiments, in which each subdirectory corresponds to figures or tables in the paper).
  


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
- Step 1. set ./experiments/Table_1_new/BrodcasTR as the working directory. 
- Step 2. run the following code to generate the synthetic data. 
```markdown
nohup  Rscript 20230918_generate_data_new_editor.R. #You will obtain the synthetic data "Simul_n1000_rep50_final_fix_new_editor.Rdata" in "./SimResults".
```
- Step 3. run the following code in the the command line.
```markdown
nohup  Rscript --vanilla "SimNonLin.R" > ./logs 2>&1 & # obtain the fitting results
nohup  Rscript --vanilla "ISE_SimNonLin.R" > ./ISE_logs 2>&1 & # compute the ISE, after you obtain the fitting results
```
<details>
  <summary>Click to view collapsible paragraph---other methods </summary>
  
#### For TLR-rescaled 
- Step 1. Put the generated synthetic data in ./experiments/Table_1_new/TLR-rescaled/SimResults
- Step 2. set ./experiments/Table_1_new/TLR-rescaled as the working directory.
- Step 3. run the following code in the the command line.
```markdown
nohup  Rscript --vanilla "SimNonLin.R" > ./logs 2>&1 &
nohup  Rscript --vanilla "ISE_SimNonLin.R" > ./ISE_logs 2>&1 &
```

#### For ENetR: 
- Step 1. Put the generated synthetic data in ./experiments/Table_1_new/ENetR /SimResults
- Step 2. set ./experiments/Table_1_new/ENetR as the working directory.
- Step 3. run the following code in the the command line.
```markdown
nohup  Rscript --vanilla "SimENetR.R" > ./logs 2>&1 &
nohup  Rscript --vanilla "ISE_SimENetR.R" > ./ISE_logs 2>&1 &
```

#### For TLR: 
- Step 1. go to ./experiments/Table_1_new/README.md and run the steps for TLR.
 </details>
 
- Remark 1: You need to install the R package snow (for the parallel computing) and the dependences of BNTR in a linux server.
- Remark 2: You can set the tuning parameters of BroadcasTR in ./experiments/Table_1_new/ParallelComput/parallel_source1000.R. for the sample size n=1000.
- Remark 3: You can set the number of CPUs for the computation in ./experiments/Table_1_new/ParallelComput/parallel_replications_big_1000K8_new.R for the sample size n=1000. 
- Remark 4: For more details, please refer to ./experiments/Table_1_new.
- Remark 5: For Fig. S.1. and Fig. S.2. in the supplementary, please refer to ./experiments/Figure_5_new and ./experiments/Figure_6_new, respectively.


## Experiment 3: The real data analysis in Section 5.2 of the main paper
For the code to reproduce Table 2 in the main paper, please refer to ./experiments/Table 2. You can obtain the results by using the following steps. 
- Step 1. Download the data from http://neurotycho.org/expdatalist/listview?task=67
- Step 2. For data preprocessing, please follow the steps in ./experiments/Table 2/ECOG_observation/README.md. You will obtain X_test_1.Rdata, X_train_1.Rdata, y_all.Rdata. 
- Step 3 Set ./experiments/Table 2 as the working directory. Create a directory named data and move X_test_1.Rdata, X_train_1.Rdata, y_all.Rdata into ./data.
- Step 4. Run "RealNonLin.R", "RealTLR-rescaled.R",and "RealENet.R" to obtain the results of BroadcasTR, TLR-rescaled, ENetR, respectively. For TLR1, please refer to ./experiments/Table 2/README.md. 

- Remark 1: You also need to install the R package snow (for the parallel computing) and the dependences of BNTR in a linux server.
- Remark 2: The data spliting follows /ECOG_observation/idtest_matrix.Rdata and ./ECOG_observation/idtrain_matrix.Rdata. 
- Remark 3: For more details, please refer to ./experiments/Table 2.



## Experiment 4: The simulated monkey electrocorticography data in Section 5.3 of the main paper





# Experiments in the paper
In summary, all the code can be found in ./experiments, in which each subdirectory corresponds to figures or tables in the paper

<details>
  <summary>Click to view collapsible paragraph---Old Version</summary>
  
## Version 3 (The previous)
### The main paper (The previous)
<details>
  <summary>Click to view collapsible paragraph---The main paper (The previous)</summary>
 
- "./experiments/Table 1":  Estimation performance in synthetic data. Reported are the averages of ISEs and the corresponding standard deviations (in parentheses) based on 50 data replications. In the first column, n is the total sample size, of which 20% were kept for the hold-out method. (X is generated from a mixture of multivariate truncated normal distribution with a Toeplitz variance matrix on the support and point mass distribution on the boundaries.)
 
- "./experiments/Table 2": Prediction performance on the monkey’s electrocorticography data. Reported are averages of MSPE and the corresponding standard deviations (in parentheses) based on 10 random splittings.
  
- "./experiments/Table 3": Prediction performance on the simulated monkey’s electrocorticography data. Reported are averages of MSPE and the corresponding standard deviations (in parentheses) based on 10 random splittings.
  
- "./experiments/Figure 4": Region selection of TLR, TLR-rescaled, ENetR, and BroadcasTR for n = 1000, of which 20% were for the hold-out method. The first column presents the true norm tensors corresponding to Cases 1–5, respectively. The rest four columns depict the estimated norm tensor with median ISEs of the comparative and proposed methods. Columns from left to right respectively correspond to TLR, TLR-rescaled, ENetR, and BroadcasTR. The plots in all columns share the same color scheme as shown in the color bar at the bottom. (X is generated from a mixture of multivariate truncated normal distribution with a Toeplitz variance matrix on the support and point mass distribution on the boundaries.)
  
- "./experiments/Figure 5": Region selection of BroadcasTR for Cases 1–5, with various sample size n = 500, 750, and 1000 (where 20% data are used for the hold-out method). All plots share the same color scheme as shown in the color bar at the bottom.(X is generated from a mixture of multivariate truncated normal distribution with a Toeplitz variance matrix on the support and point mass distribution on the boundaries.)
  
- "./experiments/Figure 6": The plot of the true and estimated entry-wise functions using BroadcasTR. From the first row to the fifth row correspond to Cases 1–5, respectively. From left to right are respectively the sample sizes n = 500, 750, and 1000. (X is generated from a mixture of multivariate truncated normal distribution with a Toeplitz variance matrix on the support and point mass distribution on the boundaries.)
  
- "./experiments/Figure 7": The estimation performance of entry-wise functions for the simulated monkey’s electrocorticography data. Each panel is the result of one comparing method. The entry is chosen to be the one with the median ISE (of BroadcasTR) among all entires, in the replicate with the median performance among 10 replicates of random splitting.
 </details>
 
### The supplementary (The previous)
<details>
  <summary>Click to view collapsible paragraph---The supplementary (The previous)</summary>
 
- "./experiments/Supp/Figure S1-S3": (1)Region selection of TLR, TLR-rescaled, ENetR, and BroadcasTR for synthetic data in Section E with n = 1000, of which 20% were for the hold-out method. The first column presents the true norm tensors corresponding to Cases 1–5, respectively. The rest four columns de- pict the estimated norm tensor with median ISEs of the comparative and proposed methods. Columns from left to right respectively correspond to TLR, TLR-rescaled, ENetR, and BroadcasTR. The plots in all columns share the same color scheme as shown in the color bar at the bottom. (2)Region selection of BroadcasTR for Cases 1–5 (from the first to the fifth row) in Section E, with various sample size n = 500, 750, and 1000 (where 20% data are used for the hold-out method). All plots share the same color scheme as shown in the color bar at the bottom. (3)The plot of the true and estimated entry-wise functions using BroadcasTR for synthetic data in Section E. From the first row to the fifth row correspond to Cases 1–5, respectively. From left to right are respectively the sample sizes n = 500, 750, and 1000. (The tensor covariate X and the error εj were generated such that Xi1,i2 ∼ Uniform[0,1])

- "./experiments/Supp/Figure S4": Estimated tensor norms of BroadcasTR for Cases 1–5 (from the first to the fifth row) in Section F where the entries of tensor covariate are Beta(2, 2) distributed. From left to right correspond to sample sizes n = 500, 750, 1000 (of which 20% data are used for the hold-out method). All plots share the same color scheme as shown in the color bar at the bottom.
  
- "./experiments/Supp/Figure S5":The LHS (LHSloss) and RHS (RHSobj) of (A.11) in Case 2 of the synthetic data in Section E. The first row show these two quantities based on 50 replications after validation, when n = 500 and n = 1000. The second row depicts the differences (RHSobj - LHSloss) accordingly.
  
- "./experiments/Supp/Figure S6-S7": (1)Prediction performance on the ADNI data. The left and right boxplots are respectively the classification accuracy of BroadcasTR and TLR-rescaled based on 10 random splittings. (2)Region selection performance on the ADNI data. The columns correspond to the slices of the tensor covariate. The rows named “Pos-” and “Neg-” are the plots of positive and negative contributions of each entry, respectively. The rows named “-TLR” and “-BroadcasTR” correspond to the tensor linear regression with rescaling strategy and the proposed broadcasted nonparametric model, respectively.
  
- "./experiments/Supp/Table S1-S3": （1）Estimation performance of synthetic data. Reported are the averages of ISEs and the corresponding standard deviations (in parentheses) based on 50 data replications. In the first column, n is the total sample size, of which 20% were kept for the hold-out method. (the tensor covariate X and the error εj were generated such that Xi1,i2 ∼ Uniform[0,1]) (2) Estimation performance of synthetic data in Section F where the entries of tensor covariate are Beta(2, 2) distributed. Reported are the averages of ISEs and the corresponding standard deviations (in parentheses) based on 50 data replications. In the first column, n is the total sample size, of which 20% were kept for the hold-out method. (3) Reported are averages of MSPE and the corresponding standard deviations (in parentheses) based on 50 replications.(GPNTE, TVGP, BroadcasTR,X = x1 ◦ x2 ∈ R64×64, where each entry of xd was independently sampled from Uniform[0,1], d = 1,2.)
 </details>

  
## Version 4 (The latest with description)
### The main paper (The latest)
- Table 1 in the main paper："./experiments/Table_1_new". Table 1. Estimation performance for the synthetic data. Reported are the averages of ISE and the corresponding standard deviations (in parentheses) based on 50 data replicates. In the first column, n is the total sample size.

- Table 2 in the main paper: "./experiments/Table 2". Table 2. Prediction performance for the monkey electrocorticography data. Reported are the averages of MSPE and the corresponding standard deviations (in parentheses) based on 10 random splits. (The monkey’s electrocorticography data is preprocessed by the file in "./experiments/Table 2/ECOG_observation". )

- Table 3 in the main paper: "./experiments/Table 3". Table 3. Prediction performance for the simulated monkey electrocorticography data. Reported are the averages of MSPE and the corresponding standard deviations (in parentheses) based on 10 random splits.

- Fig. 4 in the main paper: "./experiments/Figure_4_new". Fig. 4. Region selection of the competing methods for n = 1000. The first column presents the true norm tensors in Cases 1–5. The remaining four columns display the estimated norm tensors corresponding to the replicate of the upper median ISE performance for the competing methods. The columns from left to right correspond to TLR, TLR-rescaled, ENetR, and BroadcasTR, respectively. The plots in all columns share the same color scheme as shown in the color bar at the bottom.
  
### The supplementary (The latest)
<details>
  <summary>Click to view collapsible paragraph---The supplementary (The latest)</summary>
 
- Fig. S.1. in the supplementary: "./experiments/Figure_5_new". Fig. S.1. Region selection of BroadcasTR for the synthetic data in Section 5.1 of the main paper, with sample sizes n = 500, 750, and 1000. All plots share the same color scheme as shown in the color bar at the bottom.

- Fig. S.2. in the supplementary: "./experiments/Figure_6_new". Fig. S.2. True and estimated entry-wise functions using BroadcasTR for the synthetic data in Section 5.1 of the main paper. The five rows correspond to Cases 1–5, respectively. The columns display sample sizes n = 500, 750, and 1000.

- Fig. S.3. in the supplementary: "./experiments/Figure 7". Fig. S.3. The estimation performance of the entry-wise functions for the simulated monkey electro- corticography data in Section 5.3 of the main paper. Each panel displays one competing method. The entry is chosen to be the one with the median estimation error (of BroadcasTR) among all entries.

- Table S.1. in the supplementary: "./experiments/Table 1". Table S.1. Estimation performance for the synthetic data in Section E.2. Reported are the averages of ISE and the corresponding standard deviations (in parentheses) based on 50 data replicates. In the first column, n is the total sample size.

- Table S.2 in the supplementary: "./experiments/Supp/Table S1-S3". Table S.2. Prediction performance for the synthetic data in Section E.3. Reported are the averages of MSPE and the corresponding standard deviations (in parentheses) based on 50 replicates.

- Fig. S.4. in the supplementary: "./experiments/Figure 4". Fig. S.4. Region selection of the competing methods for the synthetic data in Section E.2 (n = 1000). The first column presents the true norm tensors in Cases 1–5. The remaining four columns display the estimated norm tensors corresponding to the replicate of the upper median ISE performance for the competing methods. The columns from left to right correspond to TLR, TLR-rescaled, ENetR, and BroadcasTR, respectively. The plots in all columns share the same color scheme as shown in the color bar at the bottom.

- Fig. S.5. in the supplementary: "./experiments/Figure 5". Fig. S.5. Region selection of BroadcasTR for Cases 1–5 (from the first to the fifth row) in Section E.2, with sample sizes n = 500, 750, and 1000. All plots share the same color scheme as shown in the color bar at the bottom.

- Fig. S.6. in the supplementary: "./experiments/Figure 6". Fig. S.6. True and estimated entry-wise functions using BroadcasTR for the synthetic data in Section E.2. The five rows correspond to Cases 1–5, respectively. The columns display sample sizes n = 500, 750, and 1000.

- Fig. S.7. in the supplementary: "./experiments/Supp/Figure S5". Fig. S.7. The LHS (LHSloss) and RHS (RHSobj) of (B.10) in Case 2 of the synthetic data in Section 5.1 of the main paper. The first row show these two quantities based on 50 replicates after vali- dation, when n = 500 and n = 1000. The second row depicts the differences (RHSobj - LHSloss) accordingly.

- Fig. S.8. in the supplementary: "./experiments/Supp/Figure S6-S7". Fig. S.8. Prediction performance for the ADNI data. The left and right boxplots are respectively the classification accuracy of BroadcasTR and TLR-rescaled based on 10 random splits.

- Fig. S.9. in the supplementary: "./experiments/Supp/Figure S6-S7". Fig. S.9. Region selection performance for the ADNI data. The columns correspond to the slices of the tensor covariate. The rows named “Pos-” and “Neg-” are the plots of positive and negative contributions of each entry, respectively. The rows named “-TLR” and “-BroadcasTR” correspond to the tensor linear regression with rescaling strategy and the proposed broadcasted nonparametric model, respectively.
 </details>

   </details>

### The main paper (The latest)
- Table 1 in the main paper："./experiments/Table_1_new".
  
- Table 2 in the main paper: "./experiments/Table 2". 

- Table 3 in the main paper: "./experiments/Table 3".
  
- Fig. 4 in the main paper: "./experiments/Figure_4_new". 
  
### The supplementary (The latest)
<details>
  <summary>Click to view collapsible paragraph---The supplementary (The latest)</summary>
 
- Fig. S.1. in the supplementary: "./experiments/Figure_5_new".
  
- Fig. S.2. in the supplementary: "./experiments/Figure_6_new". 

- Fig. S.3. in the supplementary: "./experiments/Figure 7".

- Table S.1. in the supplementary: "./experiments/Table 1".
  
- Table S.2 in the supplementary: "./experiments/Supp/Table S1-S3".
  
- Fig. S.4. in the supplementary: "./experiments/Figure 4".
  
- Fig. S.5. in the supplementary: "./experiments/Figure 5".
  
- Fig. S.6. in the supplementary: "./experiments/Figure 6".
  
- Fig. S.7. in the supplementary: "./experiments/Supp/Figure S5". 

- Fig. S.8. in the supplementary: "./experiments/Supp/Figure S6-S7". 

- Fig. S.9. in the supplementary: "./experiments/Supp/Figure S6-S7".

</details>
 

# References
- Zhou, Y., Wong, R. K. W., & He, K. (2020). Broadcasted nonparametric tensor regression. arXiv preprint arXiv:2008.12927. [\[link\]](https://arxiv.org/abs/2008.12927v2)
- Remark: there is an updated version which is under review.
