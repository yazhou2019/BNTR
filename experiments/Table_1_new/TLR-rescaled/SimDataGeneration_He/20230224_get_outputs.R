#final_fix
get_outputs <- function(v, BBprod, X_data){
  y_all = list()
  y_all[[1]]=v[[1]]+t(ctprod(BBprod[[1]],X_data,2))
  y_all[[2]]=v[[2]]+t(ctprod(BBprod[[2]],f2(X_data),2))+t(ctprod(BBprod[[8]],f4(X_data),2))+t(ctprod(BBprod[[9]],f5(X_data),2))
  y_all[[3]]=v[[3]]+t(ctprod(BBprod[[3]],f1(X_data),2))+t(ctprod(BBprod[[4]],f1(X_data),2))
  y_all[[4]]=v[[4]]+t(ctprod(BBprod[[5]],f1(X_data),2))+t(ctprod(BBprod[[6]],f3(X_data),2))
  y_all[[5]]=v[[5]]+t(ctprod(BBprod[[7]],f3(X_data),2))
  return(y_all)
}





