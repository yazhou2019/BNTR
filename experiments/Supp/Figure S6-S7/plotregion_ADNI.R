library("RNifti")
library(BNTR)

library(rTensor)
library(refund.wave)
library(refund)
library(fields)
library(stringr)



#region selection
load("./data40403.Rdata")
X_all = data$X_all
y_all = data$y_all

load("./zeroX.Rdata")
reslinear = broadcasted_sparsetenreg(X_all, y_all, r = 1, lambda = 0.075, alpha = 0.5, warmstart = 0, Replicates=5)


plus_minus_diff_linear = function(coefall, X_all){
  
  plus=0
  minus =0 
  n= dim(X_all)[4]
  for(i in 1:n){
    tem = coefall * X_all[,,,i] /2
      
      tem_plus =tem
      tem_plus[tem_plus<0] = 0
      tem_plus = array(tem_plus, dim(tem))
      tem_minus = tem
      tem_minus[tem_minus > 0]= 0
      tem_minus = array(-tem_minus, dim(tem))
      
      plus = plus + tem_plus
      minus = minus + tem_minus
  }
  res = list(plus= plus, minus =minus)
}



tem_linear = plus_minus_diff_linear(full_R(reslinear$beta),X_all)


par(mar = c(0,0,0,0),mfcol=c(3,2),mai=c(0.05,0.25,0.3,0.2))
for(slice_i in 1:3){
  po_ef_linear = tem_linear[[1]]*zeroX/774
  po_ef_linear = po_ef_linear[,,slice_i]
  po_bar_linear = seq(0,max(po_ef_linear), max(po_ef_linear)/39)
  po_bar_ef_linear = po_ef_linear
  po_bar_ef_linear[40,] = po_bar_linear
  #view(po_bar_ef)
  #image(po_bar_ef)
  library(fields)
  xx= 1:40
  yy= 1:40
  #zz = outer(xx, yy, "+")
  zz = po_ef_linear
  image.plot(xx, yy, po_ef_linear,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
  ne_ef_linear = tem_linear[[2]]*zeroX/774
  ne_ef_linear = ne_ef_linear[,,2]
  zz = ne_ef_linear
  image.plot(xx, yy, ne_ef_linear,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
}





qknots = quantile(X_all)
X_all_basis = tildePhiX_trans(X_all,knots = qknots)
res = broadcasted_sparsetenreg(X_all_basis, y_all, r = 2, lambda = 0.5, alpha = 0.5, warmstart = 1, Replicates=1)











plus_minus_diff = function(coefall, XXX,knots=stats::quantile(c(X_all), probs = c(seq(0, 1, 1/(5 - 1))))){
  
  plus=0
  minus =0 
  n= dim(XXX)[5]
  
  fhat_intten <- apply(coefall, c(1:3), fhat_int, knots = knots, order = 4)
  for(i in 1:n){
    tem = coefall *XXX[,,,,i]
    heheda = tem[,,,1]
    for(j in 2:6){
      heheda = heheda + tem[,,,j]
    }
    
    heheda = heheda 
    tem_plus = heheda
    tem_plus[tem_plus<0]= 0
    tem_plus = array(tem_plus, dim(heheda))
    tem_minus = heheda
    tem_minus[tem_minus > 0]= 0
    tem_minus = array(-tem_minus, dim(heheda))
    plus = plus + tem_plus
    minus = minus + tem_minus
  }
  plus = plus -  fhat_intten
  minus = minus - fhat_intten
  res = list(plus= plus, minus =minus)
}



  tem = full_R(res$beta)
  
  tem = plus_minus_diff(tem,X_all_basis[,,,,])

  
  par(mar = c(0,0,0,0),mfcol=c(3,2),mai=c(0.05,0.25,0.3,0.2))
  for(slice_i in 1:3){
    po_ef = tem[[1]]*zeroX/774
    po_ef = po_ef[,,slice_i]
    po_bar = seq(0,max(po_ef), max(po_ef)/39)
    po_bar_ef = po_ef
    po_bar_ef[40,] = po_bar
    library(fields)
    xx= 1:40
    yy= 1:40
    image.plot(xx, yy, po_ef,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n")
    mtext(as.character(slice_i),2,line=0.2) 
    if(slice_i==1){
      mtext("Pos-BroadcasTR",3,line=0.2) 
    }
    ne_ef = tem[[2]]*zeroX/774
    ne_ef = ne_ef[,,1]
    zz = ne_ef
    image.plot(xx, yy, ne_ef,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n")
    if(slice_i==1){
      mtext("Neg-BroadcasTR",3,line=0.2) 
    }  
  }


#plot_all


xx= 1:40
yy= 1:40
name_th = c("st", "nd", "rd")
par(mar = c(0,0,0,0),mfcol=c(3,4),mai=c(0.05,0.20,0.20,0.15))
for(slice_i in 1:3){
  po_ef_linear = tem_linear[[1]]*zeroX/774
  po_ef_linear = po_ef_linear[,,slice_i]
  image.plot(xx, yy, po_ef_linear,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
  
  mtext(str_c(slice_i,name_th[slice_i] ," slice"),2,line=0.2) 
  
  if(slice_i==1){
    mtext("Pos-TLR",3,line=0.2) 
  }
  
  ne_ef = tem[[2]]*zeroX/774
  ne_ef = ne_ef[,,1]
  zz = ne_ef
  image.plot(xx, yy, ne_ef,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n")
  if(slice_i==1){
    mtext("Neg-TLR",3,line=0.2) 
  } 
  
  po_ef = tem[[1]]*zeroX/774
  po_ef = po_ef[,,slice_i]
  image.plot(xx, yy, po_ef,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
  if(slice_i==1){
    mtext("Pos-BroadcasTR",3,line=0.2) 
  }
  
  
  ne_ef_linear = tem_linear[[2]]*zeroX/774
  ne_ef_linear = ne_ef_linear[,,2]
  image.plot(xx, yy, ne_ef_linear,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
  if(slice_i==1){
    mtext("Neg-BroadcasTR",3,line=0.2) 
  }
  
  
}




load("./data40403.Rdata")
X_tem = data$X[,,,96]

#plot_brain_alsoe
xx= 1:40
yy= 1:40
name_th = c("st", "nd", "rd")
par(mar = c(0,0,0,0),mfcol=c(5,3),mai=c(0.00,0.20,0.1,0.1))
for(slice_i in 1:3){
  
  image.plot(1:40, 1:40, X_tem[,,slice_i],col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
  mtext(str_c(slice_i,name_th[slice_i] ," slice"),3,line=0.2) 
  if(slice_i==1){
    mtext("Image",2,line=0.2) 
  }
}
  

for(slice_i in 1:3){  
  po_ef_linear = tem_linear[[1]]*zeroX/774
  po_ef_linear = po_ef_linear[,,slice_i]
  image.plot(xx, yy, po_ef_linear,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
  if(slice_i==1){
    mtext("Pos-TLR",2,line=0.2) 
  }
}

for(slice_i in 1:3){  
po_ef = tem[[1]]*zeroX/774
po_ef = po_ef[,,slice_i]
image.plot(xx, yy, po_ef,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
if(slice_i==1){
  mtext("Pos-BroadcasTR",2,line=0.2) 
}
}


for(slice_i in 1:3){ 
  ne_ef_linear = tem_linear[[2]]*zeroX/774
  ne_ef_linear = ne_ef_linear[,,2]
  image.plot(xx, yy, ne_ef_linear,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n",ann=FALSE)
  if(slice_i==1){
    mtext("Neg-TLR",2,line=0.2) 
  }
}

for(slice_i in 1:3){   
  ne_ef = tem[[2]]*zeroX/774
  ne_ef = ne_ef[,,1]
  zz = ne_ef
  image.plot(xx, yy, ne_ef,col=grey(seq(0, 1, length = 50)),xaxt='n',yaxt="n")
  if(slice_i==1){
    mtext("Neg-BroadcasTR",2,line=0.2) 
  } 
}
  
