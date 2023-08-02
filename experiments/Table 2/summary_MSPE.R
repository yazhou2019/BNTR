

# load the fitting result

val_error = c()
test_error = c()
for (i in 1:10){
  index = which.min(result[[i]][[3]][,2])
  val_error[i] = result[[i]][[3]][index,2]
  test_error[i] = result[[i]][[3]][index,3]
}

