
#### Create lag matrix of X ####

library(midasr)

# input an object of class ts

create_lag_Xmatrix <- function(tsdata,
                               lags,
                               frequency){
  X_matrix <- mls(tsdata, lags, frequency)
  compile_names <- c()
  for(i in 0:(ncol(X_matrix)-1)){
    compile_names <- c(compile_names,paste0("lag",i))
  }
  colnames(X_matrix) <- compile_names
  
  return(X_matrix)
}



fit_Minla <- function(){
  
  
  
}
