# Data
- The real data is preprocessed by the file in ./ECOG_observation.
- The sfSource() of (BroadcasTR, TLR-rescaled, ENetR) will use a R file including loading (X_test_1.Rdata, X_train_1.Rdata, y_all.Rdata, idtest_matrix.Rdata, idtrain_matrix.Rdata). The data can be put in the file named data.

# TLR1
- The code depends on the MATLAB code in https://hua-zhou.github.io/TensorReg/
- The real data in .Rdata format need to be transformed to the .mat format; see ./TLR1/Data_generation_to_MATLAB_ECoG.R.
- X_train_1.Rdata, X_test_1.Rdata, y_all.Rdata need to put in the file ./TLR1/data


# BroadcasTR, TLR-rescaled, ENetR
- See RealNonLin.R, RealTLR-rescaled.R,and RealENet.R, respectively


# MSPE
- The MSPE of BroadcasTR, TLR-rescaled, ENetR can be obtain by loading the fitting results (.Rdata) and utilizing the code in summary_MSPE.R
- The MSPE of TLR can be obtain from the fitting results. 


