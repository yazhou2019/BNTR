# Data
- The real data is preprocessed by the file in ./generate_real_sim_Data_new.R.
- Before run generate_real_sim_Data_new.R, you need to put X_test_1.Rdata, X_test_1.Rdata, y_all.Rdata and monkey_1_new_K11_tuning.Rdata in ./data. These .Rdata are the same as those of Table 2


# TLR1
- The code depends on the MATLAB code in https://hua-zhou.github.io/TensorReg/
- - X_train_1.Rdata, X_test_1.Rdata, y_all_sim_new2.Rdata need to put in the file ./TLR1/data
- The real data in .Rdata format need to be transformed to the .mat format by ./TLR1/Data_generation_to_MATLAB_ECoGSim.R
- run TLR_for_RealDataSim.m, and you will obtain the fitting result.


# BroadcasTR, TLR-rescaled, ENetR
- See RealNonLin_repro_K11_tuning.R, RealTLR-rescaled.R, and RealENet.R, respectively

# MSPE
- The MSPE of BroadcasTR, TLR-rescaled, ENetR can be obtain by loading the fitting results (.Rdata) and utilizing the code in summary_MSPE.R
- The MSPE of TLR can be obtain from the fitting results. 


 


