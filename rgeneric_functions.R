

rgeneric.Beta.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                       "log.prior", "quit"),
                               theta = NULL){
  envir = parent.env(environment())
  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)
  
  interpret.theta = function() {
    
    lag_k <- ncol(x)-2
    gamma2 <- exp(theta[2L]) + 1
    for(lag in 0:lag_k){
      x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      temp <- 1 + (1-x_temp)^(gamma2-1)
      assign(paste0("psi",lag), temp)
      
    }
    
    compile_sum <- 0
    for(lag in 0:lag_k){
      compile_sum <- compile_sum + get(paste0("psi", lag))
    }
    
    for(lag in 0:lag_k){
      temp <- get(paste0("psi",lag))
      assign(paste0("w",lag),temp/compile_sum)
    }
    
    out_list <- vector(mode = "list", length = 1 + lag_k + 1)
    out_list[[1]] <- theta[1L]
    for(lag in 0:lag_k){
      out_list[[lag+2]] <- get(paste0("w",lag))
    }
    
    compile_names <- c()
    for(lag in 0:lag_k){
      compile_names <- c(compile_names, paste0("w",lag))
    }
    names_vec <- c("beta1", compile_names)
    names(out_list) <- names_vec
    return(out_list)
  }
  graph = function() {
    G = Diagonal(n = length(x$lag0), x=1) 
    return(G)
  }
  Q = function() {
    Q = prec.high * graph() 
    return(Q)
  }
  mu = function() {
    par = interpret.theta() 
    
    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }
    
    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }
    
    return(par$beta1 * agg)
  }
  log.norm.const = function() { 
    return(numeric(0))
  }
  log.prior = function() {
    par = interpret.theta()
    n.params <- length(par)
    val <- 0
    for(i in 1:n.params){
      temp <- dnorm(par[[i]], mean=0, sd=1, log=TRUE)
      val <- val + temp
    }
  }
  initial = function() { 
    return(rep(1, 2))
  }
  quit = function() { 
    return(invisible())
  }
  
  val = do.call(match.arg(cmd), args = list())
  return(val)
}



rgeneric.Almon.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                        "log.prior", "quit"),
                                theta = NULL){
  envir = parent.env(environment())
  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)
  
  interpret.theta = function() {
    
    lag_k <- ncol(x)-2
    gamma1 <- 0.01*sin(theta[2L])
    gamma2 <- 0.01*sin(theta[3L])
    for(lag in 0:lag_k){
      temp <- exp(gamma1*(lag^1) + gamma2*(lag^2))
      assign(paste0("psi",lag), temp)
    }
    
    compile_sum <- 0
    for(lag in 0:lag_k){
      compile_sum <- compile_sum + get(paste0("psi", lag))
    }
    
    for(lag in 0:lag_k){
      temp <- get(paste0("psi",lag))
      assign(paste0("w",lag),temp/compile_sum)
    }
    
    out_list <- vector(mode = "list", length = 1 + lag_k + 1)
    out_list[[1]] <- theta[1L]
    for(lag in 0:lag_k){
      out_list[[lag+2]] <- get(paste0("w",lag))
    }
    
    compile_names <- c()
    for(lag in 0:lag_k){
      compile_names <- c(compile_names, paste0("w",lag))
    }
    names_vec <- c("beta1", compile_names)
    names(out_list) <- names_vec
    return(out_list)
    
  }
  graph = function() {
    G = Diagonal(n = length(x$lag0), x=1) 
    return(G)
  }
  Q = function() {
    Q = prec.high * graph() 
    return(Q)
  }
  mu = function() {
    par = interpret.theta() 
    
    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }
    
    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }
    
    return(par$beta1 * agg)
  }
  log.norm.const = function() { 
    return(numeric(0))
  }
  log.prior = function() {
    par = interpret.theta()
    n.params <- length(par)
    val <- 0
    for(i in 1:n.params){
      temp <- dnorm(par[[i]], mean=0, sd=1, log=TRUE)
      val <- val + temp
    }
  }
  initial = function() { 
    return(rep(1, 3))
  }
  quit = function() { 
    return(invisible())
  }
  
  val = do.call(match.arg(cmd), args = list())
  return(val)
}



rgeneric.Hyperbolic.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                             "log.prior", "quit"),
                                     theta = NULL){
  envir = parent.env(environment())
  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)
  
  interpret.theta = function() {
    
    lag_k = ncol(x)-2
    gamma_val = exp(theta[2L])/(1+exp(theta[2L]))
    for(lag in 0:lag_k){
      temp <- gamma(lag+gamma_val) / (gamma(lag+1)*gamma(gamma_val))
      assign(paste0("psi",lag), temp)
    }
    
    compile_sum <- 0
    for(lag in 0:lag_k){
      compile_sum <- compile_sum + get(paste0("psi", lag))
    }
    
    for(lag in 0:lag_k){
      temp <- get(paste0("psi",lag))
      assign(paste0("w",lag),temp/compile_sum)
    }
    
    out_list <- vector(mode = "list", length = 1 + lag_k + 1)
    out_list[[1]] <- theta[1L]
    for(lag in 0:lag_k){
      out_list[[lag+2]] <- get(paste0("w",lag))
    }
    
    compile_names <- c()
    for(lag in 0:lag_k){
      compile_names <- c(compile_names, paste0("w",lag))
    }
    names_vec <- c("beta1", compile_names)
    names(out_list) <- names_vec
    return(out_list)
    
  }
  graph = function() {
    G = Diagonal(n = length(x$lag0), x=1) 
    return(G)
  }
  Q = function() {
    Q = prec.high * graph() 
    return(Q)
  }
  mu = function() {
    par = interpret.theta() 
    
    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }
    
    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }
    
    return(par$beta1 * agg)
  }
  log.norm.const = function() { 
    return(numeric(0))
  }
  log.prior = function() {
    par = interpret.theta()
    n.params <- length(par)
    val <- 0
    for(i in 1:n.params){
      temp <- dnorm(par[[i]], mean=0, sd=1, log=TRUE)
      val <- val + temp
    }
  }
  initial = function() { 
    return(rep(1, 2))
  }
  quit = function() { 
    return(invisible())
  }
  10
  val = do.call(match.arg(cmd), args = list())
  return(val)
}
