# Data
- You need to put the fitting results and ise of BroadcasTR in SimResults (including the case n=500, 750, 1000). Those .Rdata files can be obtain by using the codes of Table_1_new. 


# Region selection of BroadcasTR for Cases 1–5, with various sample size n = 500, 750, and 1000 
- Run 20230918_plot_vary.R, you will get the plot results; the dependence can be founded in "/experiments/Table_1_new". You can revise the code in 20230918_plot_vary.R as follow, 
replace
```
setwd("~/Desktop/Research/JRSSB/换editor上传/Table_1/BroadcasTR")
```
by
```
setwd("./Table_1_new/BroadcasTR")
```
and replace
```
setwd("~/Desktop/Research/JRSSB/换editor上传/Table_1")
```
by
```
setwd("./Table_1_new")
```
