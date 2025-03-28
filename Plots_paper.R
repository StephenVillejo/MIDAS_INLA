#### Plots of weights for 3 constraint functions ####

##### Almon #####

vals <- vector(mode = "list", length = 2)
vals[[1]]$gamma1 = 0.002 # rapidly declining weights
vals[[1]]$gamma2 = -.005
vals[[2]]$gamma1 = 0.001 # slowly declining weights
vals[[2]]$gamma2 = -.0005

compile_weights <- vector(mode = "list", length = 2)
for(j in 1:2){
  for(i in 0:100){
    compile_weights[[j]] <- c(compile_weights[[j]],
                              exp((vals[[j]]$gamma1*(i^1)) + 
                                    vals[[j]]$gamma2*(i^2))) 
  }
  compile_weights[[j]] <- compile_weights[[j]]/sum(compile_weights[[j]])
}



par(mfrow=c(1,2))
plot(compile_weights[[1]], type="l", col="blue", 
     ylab="weights", xlab="lag", 
     main=expression(paste(gamma[1]==0.002, " ",gamma[1]==-0.005)))
plot(compile_weights[[2]], type="l", col="blue", 
     ylab="weights", xlab="lag", 
     main=expression(paste(gamma[1]==0.001, " ",gamma[1]==-0.0005)))




##### Beta #####

gamma1 = 1
gamma2 = 4

gamma1 = 1
gamma2 = 8

compile_weights <- c()
for(i in 0:100){
  x <- 0.0001 + (1-0.0001)*((i-1)/(100-1))
  compile_weights <- c(compile_weights,
                       (x^(gamma1-1))*((1-x))^(gamma2-1)) 
}
compile_weights <- compile_weights/sum(compile_weights)
plot(compile_weights)


##### Hyperbolic scheme #####

gamma <- .5

gamma <- .9

compile_weights <- c()
for(lag in 0:100){
  compile_weights <- c(compile_weights,
                       gamma(lag+gamma) / (gamma(lag+1)*gamma(gamma)))
}
compile_weights <- compile_weights/sum(compile_weights)
plot(compile_weights)




