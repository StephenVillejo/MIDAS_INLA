source("rgeneric_functions.R")
source("estimation_functions.R")

library(midasr)
library(ggplot2)
library(INLA)


#### Simulated data ####

n <- 600 # 50 years of monthly data
beta0 = 10
trend = .05
beta1 = 3
sigma_e = 1

gamma_param <- 0.3
lag_k <- 60
xi <- c()
for(lag in 0:lag_k){
  temp <- gamma(lag+gamma_param)/(gamma(lag+1)*gamma(gamma_param))
  xi <- c(xi, temp)
}
weights <- xi/sum(xi)
plot(weights, type = "l")

set.seed(120054)
x <- rnorm(30*n) 
y <- beta0 + trend*c(1:n) + beta1*mls(x,0:lag_k,30)%*%weights + rnorm(n, sd = sigma_e)
plot(y, type = "l")
plot(x, typ = "l")


y_ts <- ts(y, frequency = 1)
x_ts <- ts(x, frequency = 30)
png("/Users/stephenjunvillejo/Downloads/SimHighFrequency_hyperbolic.png", width=17, height=10, units = 'cm', res = 300)
autoplot(x_ts, color = "brown") + theme_bw() + geom_point(color = "red", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()
png("/Users/stephenjunvillejo/Downloads/SimLowFrequency_hyperbolic.png", width=17, height=10, units = 'cm', res = 300)
autoplot(y_ts, color = "black") + theme_bw() + geom_point(color = "gray", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()



#### Model fitting ####

Midas_objects <- fit_Minla(xdata = x_ts,
                           ydata = y_ts,
                           constraint = "hyperbolic",
                           K = 0:60,
                           m = 30,
                           lagY = 0) 

res = inla(y ~ 1 + f(idx, model = Midas_objects$rgen, n = nrow(Midas_objects$data)) +
             f(trend, model = "linear"),
           data = data.frame(y = Midas_objects$data$y, 
                             trend = 1:nrow(Midas_objects$data),
                             idx = 1:nrow(Midas_objects$data)),
           verbose = TRUE,
           control.family = list(hyper = list(prec = list(prior = "pc.prec", param = c(sigma_e, 0.5)))),
           control.compute=list(config = TRUE))

summary(res)


#### Compare observed versus predicted values ####

pred_data <- Midas_objects$data
pred_res <- predict_midas(model = res,
                          data = pred_data)
both_ts <- ts(data.frame(observed = as.vector(Midas_objects$data$y),
                         predicted = pred_res$preds$mean),
              start = 3, end = 598)

png("/Users/stephenjunvillejo/Downloads/PredsVsObs_hyperbolic.png", width=25, height=12, units = 'cm', res = 300)
autoplot(both_ts) +
  theme_bw() +
  scale_colour_manual(
    values = c("blue","red")
  ) +
  ylab("y") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()



#### Plot posterior marginals ####

png("/Users/stephenjunvillejo/Downloads/ParamEstimates_hyperbolic.png", width=30, height=10, units = 'cm', res = 300)
par(mfrow=c(1,4))

plot(inla.smarginal(res$marginals.fixed[["(Intercept)"]]),
     type="l", lwd=3, col="red", xlab=expression(beta[0]), ylab="",
     cex.lab = 1.9, cex.axis=1.5)
abline(v = beta0, col = 'blue', lty=1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed$`(Intercept)`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed$`(Intercept)`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Theta1 for idx']]),
     type="l", lwd=3, col="red", xlab=expression(beta[1]), ylab="",
     cex.lab = 1.9, cex.axis=1.5)
abline(v = beta1, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta1 for idx`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Theta1 for idx`), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.fixed[['trend']]),
     type="l", lwd=3, col="red", xlab=expression(beta[2]), ylab="",
     cex.lab = 1.9, cex.axis=1.5)
abline(v = trend, col = 'blue', lty = 1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed[['trend']]), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.fixed[['trend']]), prob = 0.975), lty = 2)

plot(inla.smarginal(res$marginals.hyperpar[['Precision for the Gaussian observations']]),
     type="l", lwd=3, col="red", xlab=expression(1/sigma[epsilon]^2), ylab="",
     cex.lab = 1.5, cex.axis=1.5)
abline(v = 1/(sigma_e^2), col = 'blue', lty=1, lwd = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Precision for the Gaussian observations`), prob = 0.025), lty = 2)
abline(v = quantile(inla.rmarginal(200, res$marginals.hyperpar$`Precision for the Gaussian observations`), prob = 0.975), lty = 2)

dev.off()


#### Compare true and estimated (lag) weights ####

theta2_samples <- inla.rmarginal(200, res$marginals.hyperpar$`Theta2 for idx`)
gamma_samples <- exp(theta2_samples)/(1+exp(theta2_samples))
compile_sum <- vector(length = 200)
for(lag in 0:lag_k){
  temp <- gamma(lag+gamma_samples)/(gamma(lag+1)*gamma(gamma_samples))
  assign(paste0("psi_sample_",lag),temp)
  compile_sum <- compile_sum + temp
}
for(lag in 0:lag_k){
  temp <- get(paste0("psi_sample_", lag))
  assign(paste0("w",lag),temp/compile_sum)
}
compile_w_list <- vector("list", length = lag_k + 1)
for(i in 0:lag_k){
  temp <- get(paste0("w",i))
  compile_w_list[[i+1]] <- temp
}
mean <- lapply(1:(lag_k+1), function(x) mean(compile_w_list[[x]]))
q2.5 <- lapply(1:(lag_k+1), function(x) quantile(compile_w_list[[x]],probs=0.025))
q97.5 <- lapply(1:(lag_k+1), function(x) quantile(compile_w_list[[x]],probs=0.975))

weights <- weights[(lag_k+1):1]
weights_ests_df <- data.frame(true = weights,
                              mean = unlist(mean),
                              q2.5 = unlist(q2.5),
                              q97.5 = unlist(q97.5),
                              lag = 0:lag_k)

weights_ests_df$lag <- factor(weights_ests_df$lag, levels=unique(weights_ests_df$lag))


png("/Users/stephenjunvillejo/Downloads/WeightsEstimates_hyperbolic.png", width=35, height=15, units = 'cm', res = 300)
ggplot(weights_ests_df, aes(x=lag, y=mean)) + 
  #geom_function(fun = weights_func, xlim = c(1, 40), color = "gray", n = 500) +
  geom_point(aes(col="Posterior mean")) +
  geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=.2,
                position=position_dodge(0.05), col = "black") + 
  geom_point(aes(x=lag, y=true, color = "True value")) + 
  scale_color_manual(name='',
                     breaks=c('Posterior mean', 'True value'),
                     values=c('Posterior mean'='red', 'True value'='blue')) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank()) + 
  ylab("")
dev.off()




