
#### Create lag matrix of X ####

library(midasr)
library(rlist)
library(matrixStats)

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



fit_Minla <- function(xdata,
                      ydata,
                      constraint,
                      K,
                      m,
                      lagY = 0){
  
  X_matrix <- create_lag_Xmatrix(tsdata = xdata,
                                 lags = K,
                                 frequency = m)
  
  temp_data <- as.data.frame(X_matrix)
  
  rm.row <- which(complete.cases(temp_data) == FALSE)
  if(length(rm.row > 0) > 0){
    temp_data$y <- as.vector(ydata)
    temp_data <- temp_data[-rm.row,]
  }else{
    temp_data$y <- as.vector(ydata)
  }

  
  if(constraint == "beta"){
    rgen = inla.rgeneric.define(model = rgeneric.Beta.midas,
                                x = temp_data)
  }else if(constraint == "beta2"){
    rgen = inla.rgeneric.define(model = rgeneric.Beta2.midas,
                                x = temp_data)
  }else if(constraint == "almon2"){
    rgen = inla.rgeneric.define(model = rgeneric.Almon2.midas,
                                x = temp_data)
  }else if(constraint == "almon3"){
    rgen = inla.rgeneric.define(model = rgeneric.Almon3.midas,
                                x = temp_data)
  }else if(constraint == "hyperbolic"){
    rgen = inla.rgeneric.define(model = rgeneric.Hyperbolic.midas,
                                x = temp_data)
  }
  
  if(lagY == 0){
    data = temp_data
  }else{
    lagYdata <- matrix(NA, nrow = nrow(temp_data), ncol = lagY)
    for(i in 1:lagY){
      lagYdata[,i] <- as.vector(mls(temp_data$y, i, 1))
    }
    lagYdata <- as.data.frame(lagYdata)
    compile_names <- c()
    for(i in 1:lagY){
      compile_names <- c(compile_names,paste0("lagy",i))
    }
    names(lagYdata) <- compile_names
    data <- cbind(temp_data, lagYdata)
  }
  
  return(out = list(data = data,
                    X_matrix = X_matrix,
                    rgen = rgen,
                    rm.row = rm.row))
  
}


predict_midas <- function(model,
                          data){
  
  compile_preds <- lapply(1:10, function(x){
    temp <- inla.posterior.sample(n = 1, res_Beta)
    return(temp[[1]]$latent[1:nrow(data)])
  })
  compile_preds <- list.cbind(compile_preds)
  computed_preds <- list(mean = rowMeans(compile_preds),
                         sd = apply(compile_preds, 1, sd),
                         q2.5 = rowQuantiles(compile_preds, probs = 0.025),
                         q97.5 = rowQuantiles(compile_preds, probs = 0.975))
  
  return(out = list(preds = computed_preds))
  
}

                            
