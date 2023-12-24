#!/bin/bash

#SBATCH -J Nonlinear_sim

### 队列名，不变
#SBATCH -p amd_512

### 节点数，不变
#SBATCH -N 1

### 线程数，可以理解为核心数
### 可改，一般16，最大32
#SBATCH -n 17

#SBATCH --mem-per-cpu=12G

######SBATCH  --time=12:00:00


###conda activate R-420

Rscript SimNonLin.R

##nohup Rscript --vanilla SimNonLin.R  > /public3/home/scg5453/yzhou/BNTR/FullSimPaper/simul_BNTR_1000_K8 2>&1 &
