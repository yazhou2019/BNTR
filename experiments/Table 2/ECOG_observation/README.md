# Data preprocess for ECOG data.
- Step 1. Download the data from http://neurotycho.org/expdatalist/listview?task=67
- Step 2. Run "./DataPrepro/01Read_preprocessing.R" in R
- Step 3. Run "./DataPrepro/02Wavenet_preprocessing.m" in MATLAB
- Step 4. Run "./DataPrepro/03aTensorize_preprocessing.R" and "03boutput_preprocessing.R" in R
You will get "X_train_1.Rdata", "X_test_1.Rdata" and "y_all.Rdata"


