# Simulation setup

if(!interactive()){
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
}else{
  args <- c(100)
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# install.packages( "MatrixModels", type="win.binary" )  

library(nleqslv)
library(CVXR)
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
library(caret)
library(tidyverse)
library(kableExtra)
library(truncnorm)  # For dtruncnorm and normalization

set.seed(11)
SIMNUM = args[1]

if(!interactive()){
  dir.create(timenow0)
  setwd(timenow0)
  
  sink(timenow, append=TRUE)
}
# setup parallel backend to use many processors
cores = min(detectCores() - 3, 101)
print(paste("cores =", cores))
# cl <- makeCluster(cores, outfile = timenow) #not to overload your computer
cl <- makeCluster(cores)
registerDoParallel(cl)
registerDoRNG(seed = 11)

# modelcases = expand.grid(c(T,F),c(T,F),c(T,F))
modelcases = expand.grid(c(T,F),c(T,F))
# modelcases = expand.grid(c(T,F))
rnames <- apply(apply(modelcases, 2, ifelse, "C", "M"), 1, paste, collapse = "")

# theta = 210 # Kang & Schaffer

final_res <- foreach(
  simnum = 1:SIMNUM, 
  .packages = c("nleqslv", "CVXR", "caret", "truncnorm"),
  .errorhandling="pass") %dopar% {
    theta_mat = NULL
    var_mat = NULL
    CR_mat = NULL
    alpha_vec = NULL
    for(cnt in 1:nrow(modelcases)){
      # RM = modelcases[cnt, 1]
      OM = modelcases[cnt, 2]
      # VM = F
      VM = modelcases[cnt, 1]
      
      # Simulation setup of
      # # Kang & Schafer, 2007
      # n = 1000; # n = 5000
      # pcol = 4; # pcol = 4 is True
      # x = matrix(rnorm(n * pcol, 0, 1), nc= pcol)
      # z = cbind(exp(x[,1] /  2), x[,2] / (1 + exp(x[,1])) + 10,
      #           (x[,1] * x[,3] / 25 + 0.6)^3, (x[,2] + x[,4] + 20)^2) # True
      # # z = cbind(exp(x[,1] /  2), x[,2] / (1 + exp(x[,1])) + 10,
      # #           (x[,3] * x[,4] / 25 + 0.6)^3, (x[,3] + x[,4] + 20)^2)
      # pi = 1 / (1 + exp(-(-x[,1] + 0.5 * x[,2] -
      #                       0.25 * x[,3] - 0.1 * x[,4]))) # True
      # # pi = 1 / (1 + exp(-(-x[,1] + 0.5 * x[,2]))) 
      # vx = rep(160, n)
      # # vx = rep(1, n) # True
      # # VM = T
      # # if(VM == T){
      # #   vx = rep(160, n)
      # # }else{
      # # vx = (x[,1]^2 + exp(x[,3]) / 1.65) * 80
      # #   # vx = ifelse(x[,1]^2 > 2, 2, x[,1]^2)
      # # }
      # e = rnorm(n, 0, sqrt(vx))
      # y = theta + 27.4 * x[,1] + 13.7 * x[,2] +
      #   13.7 * x[,3]+ 13.7 * x[,4] + e # True
      # # y = theta + 13.7 * x[,3] + 13.7 * x[,4] + e
      
      # # Simulation setup of
      # # Wang & Kim 2024
      # n = 2000 # n = 200
      # pcol = 4; # pcol = 4 is True
      # x = matrix(rnorm(n * pcol, 2, 1), nc= pcol)
      # z = x
      # if(RM  == T){
      #   pi = 1 / (1 + exp(-(1 - x[,1] + 0.5 * x[,2] + 0.5 * x[,3] - 0.25 * x[,4])))
      # }else if(RM == F){
      #   pi = 1 / (1 + exp(-(1 - 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
      #                         0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
      # }
      # vx = rep(1, n) # original
      # e = rnorm(n, 0, sqrt(vx))
      # 
      # if(OM == T){
      #   y = 1 + x[,1] + x[,2] + x[,3] + x[,4] + e
      #   theta = 9
      # }else if(OM == F){
      #   y = 1 + 0.5 * x[,1]^2 * x[,2]^2 + 0.2 * x[,2]^2 * x[,3]^3 +
      #     0.5 * cos(x[,3]^2) + 0.5 * sin(x[,4]^2) + e
      #   theta = 27.58047  # To be corrected
      # }
      
      # # Simulation setup of
      # # Wang & Kim 2024, modified
      # n = 2000 # n = 200
      # pcol = 4; # pcol = 4 is True
      # x = matrix(rnorm(n * pcol, 2, 1), nc= pcol)
      # z = x
      # if(RM  == T){
      #   pi = 1 / (1 + exp(-(1 - x[,1] + 0.5 * x[,2] + 0.5 * x[,3] - 0.25 * x[,4])))
      # }else if(RM == F){
      #   pi = 1 / (1 + exp(-(1 - 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
      #                         0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
      # }
      # vx = rep(1, n) # original
      # e = rnorm(n, 0, sqrt(vx))
      # 
      # if(OM == T){
      #   y = 1 + x[,1] - x[,2] + x[,3] - x[,4] + e
      #   theta = 1
      # }else if(OM == F){
      #   y = 1 - 0.5 * x[,1]^2 * x[,2]^2 + 0.2 * x[,2]^2 * x[,3]^3 -
      #     0.5 * cos(x[,3]^2) + 0.5 * sin(x[,4]^2) + e
      #   theta = 2.551373  # To be corrected
      # }
      
      # # Simulation setup of
      # # Wang & Kim 2024, modified
      # n = 2000 # n = 200
      # pcol = 4; # pcol = 4 is True
      # x = matrix(rnorm(n * pcol, 2, 1), nc= pcol)
      # # x = ifelse(x < 0, 0, x)
      # # x = ifelse(x > 4, 4, x)
      # z = x
      # if(RM  == T){
      #   pi = 1 / (1 + exp(-(-0.5 - x[,1] + 0.5 * x[,2] + 0.5 * x[,3] - 0.25 * x[,4])))
      # }else if(RM == F){
      #   #   pi = 1 / (1 + exp(-(- 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.75 * (x[,2] - 2)^2 -
      #   #                         0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
      #   # pi = ifelse(pi < 0.05, 0.05, pi)
      #   # pi = ifelse(pi > 0.7, 0.7, pi)
      #     # pi = 1 / (1 + exp(-(1 - 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
      #     #                       0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
      #     # pi = pt(- 0.1 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
      #     #           0.2 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2, 3)
      #   pi = 1 / (1 + exp(-(- 0.1 * x[,2] * log(abs(x[,1] + 5)) - 0.1 * exp(x[,3] - 2) -
      #                         0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
      #     pi = ifelse(pi < 0.05, 0.05, pi)
      #     # pi = ifelse(pi > 0.7, 0.7, pi)
      #     summary(pi)
      #     mean(pi)
      # }
      # vx = rep(2, n) # original
      # e = rnorm(n, 0, sqrt(vx))
      # 
      # if(OM == T){
      #   y = 1 + x[,1] - x[,2] + x[,3] - x[,4] + e
      #   theta = 1
      # }else if(OM == F){
      #   # y = 1 - 0.5 * (x[,1] - 2)^2 * (x[,2]-5)^2 + 0.2 * (x[,2] - 1)^2 * (x[,3] - 1.5)^3 -
      #   #   0.5 * cos(x[,3]^2) * sin(x[,4]^2) + e
      #   # theta = -3.351918
      #     # y = 1 + 0.5 * x[,1]^2 * x[,2]^2 + 0.2 * x[,2]^2 * x[,3]^3 +
      #     #   0.5 * cos(x[,3]^2) + 0.5 * sin(x[,4]^2) + e
      #     # theta = 27.58047
      #     # y = 1 + 0.5 * x[,1]^2 * x[,2]^2 + 0.2 * x[,2]^2 * x[,3]^2 + 0.5 * sin(x[,4]^2) + e
      #     # theta = 18.5629
      #     
      #     # y = 1 - 0.2 * x[,1]^4 + 0.2 * x[,3]^4 + e
      #     # theta = 1
      #     
      #     # y = 1 - 0.2 * (x[,1] - 2)^4 + 0.2 * (x[,3] - 2)^4 + e
      #     # theta = 1
      # 
      #   # y = 1 - 0.5 * (x[,1] - 2)^2 + 0.2 * (x[,2] - 1)^2 + 0.2 * (x[,3] - 1.5)^2 - 0.2 * (x[,4] - 2.5)^2  + e
      #   # theta = 0.9  # OM1
      #   
      #   # y = 1 - 0.5 * (x[,1] - 2)^2 + 0.2 * (x[,2] - 1)^2 * (x[,3] - 1.5)^3 -
      #   #     0.5 * cos(x[,3]^2) * sin(x[,4]^2)  + e
      #   # theta = 1.148082 
      # 
      #   # y = 1 - 0.5 * (x[,1] - 3)^2 * (x[,2]-5)^2 + 0.2 * (x[,2] - 1)^2 * (x[,3] - 1.5)^3 -
      #   #   0.5 * cos(x[,3]^2) - 0.5 * sin(x[,4]^2) + e
      #   # theta = -8.430471
      # 
      #   y = 1 - 0.5 * (x[,1] - 2)^2 * (x[,2]-5)^2 + 0.2 * (x[,2] - 1)^2 * (x[,3] - 1.5)^3 -
      #     0.5 * cos(x[,3]^2) * sin(x[,4]^2) + e # OM2
      #   theta = -3.351918
      #   
      #   # y = 1 - 0.1 * (x[,1] - 2)^2 * (x[,2]-5)^2 + 0.02 * (x[,2] - 1)^2 * (x[,3] - 1.5)^3 -
      #   #   0.5 * cos(x[,3]^2) * sin(x[,4]^2) + e
      #   # theta = 0.06308178 
      #   # boxplot(y)
      #   # hist(y)
      # }
      
      # # Simulation setup of
      # # Wang & Kim 2024, modified
      # # n = 3000 # n = 200
      # n = 4000 # n = 200
      # pcol = 6; # pcol = 4 is True
      # x = matrix(rnorm(n * pcol, 2, 1), nc= pcol)
      # # x = matrix(rtruncnorm(n * pcol, a=0, b=4, mean=2, sd=1), ncol=pcol)
      # # var(x)
      # z = x
      # if(RM  == T){
      #   pi = 1 / (1 + exp(-(-0.5 - x[,3] + 0.5 * x[,4] + 0.5 * x[,5] - 0.25 * x[,6])))
      # }else if(RM == F){
      #   #   pi = 1 / (1 + exp(-(- 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.75 * (x[,2] - 2)^2 -
      #   #                         0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
      #   # pi = ifelse(pi < 0.05, 0.05, pi)
      #   # pi = ifelse(pi > 0.7, 0.7, pi)
      #   pi = 1 / (1 + exp(-(-0.5 - 0.5 * (x[,3] - 2) * (x[,4] - 5) - 0.5 * (x[,4] - 2.5)^2 -
      #                         0.5 * (x[,5] - 1.5)^2 * (x[,6] - 2)^2 )))
      #   # pi = pt(- 0.1 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
      #   #           0.2 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2, 3)
      #   # pi = 1 / (1 + exp(-(- 0.1 * x[,4] * log(abs(x[,3] + 5)) - 0.1 * exp(x[,5] - 2) -
      #   #                       0.5 * (x[,5] - 1.5)^2 * (x[,6] - 2)^2 )))
      #   pi = ifelse(pi < 0.05, 0.05, pi)
      #   # pi = ifelse(pi > 0.7, 0.7, pi)
      #   summary(pi)
      #   mean(pi)
      # }
      # vx = rep(1, n) # original
      # # e = rnorm(n, 0, 1)
      # e = rnorm(n, 0, sqrt(apply(cbind(0.5, x[,c(2,3)]^2/ 4), 1, max)))
      # 
      # if(OM == T){
      #   y = 1 + x[,1] - x[,2] + x[,3] - x[,4] + e
      #   theta = 1
      # }else if(OM == F){
      #   # y = 1 - 0.5 * (x[,1] - 2)^2 + 0.2 * (x[,2] - 1)^2 + 0.2 * (x[,3] - 1.5)^2 - 0.2 * (x[,4] - 2.5)^2  + e
      #   # theta = 0.9  # OM1
      #   
      #   y = 1 - (x[,1] - 2) * (x[,2] - 1) + 0.75 * (x[,3] - 1.5)^2 - sin(x[,4]^2)  + e
      #   theta = 1.805657  # OM1
      #   
      #   # q1y = quantile(y, .25)
      #   # q3y = quantile(y, .75)
      # 
      #   # y = ifelse(y < -2, -2, y)
      #   # y = ifelse(y > 4, 4, y)
      #   
      #   # y = 1 - 0.1 * (x[,1] - 2)^2 * (x[,2]-5)^2 + 0.05 * (x[,3] - 1)^2 * (x[,4] - 1.5)^3 -
      #   #   0.5 * cos(x[,3]^2) * sin(x[,4]^2) + e # OM2
      #   # theta = 0.1605818 
      #   # boxplot(y, horizontal = T)
      # }
      # # theta = mean(y)
      
      # Sim ogn
      n    <- 1000
      pcol <- 3
      x    <- matrix(rnorm(n * pcol, mean = 2, sd = 1), ncol = pcol)
      # colnames(x) <- c("X1","X2")
      z    <- x    # no distinction

      # c  <- -0.5 + 0.3 * (x[,2] - 2)^2 - 0.2 * (x[,3] - 2)^2
      # pi <- 1 / (1 + exp(-c))
      
      # ## Response model (RM)
      # if (RM) {
      #   # Correct RM: simple logistic in X1, X2
      #   c  <- -0.5 - 0.25 * x[,2] + 0.5 * x[,3]
      #   pi <- 1 / (1 + exp(-c))
      # } else {
      #   # Misspecified RM: nonlinear in X1, X2
      #   c  <- -0.5 + 0.3 * (x[,2] - 2)^2 - 0.2 * (x[,3] - 2)^2
      #   pi <- 1 / (1 + exp(-c))
      # }
      # delta = rbinom(n, 1, pi)

      ## Error term
      vx = rep(1, n) # original
      if(!VM){
        e = rnorm(n, 0, sqrt(apply(cbind(0.5, x[,c(2,3)]^2/ 4), 1, max)))
      }else{
        e <- rnorm(n, 0, 1)        
      }

      ## Outcome model (OM)
      if (OM) {
        # Correct OM: linear in X1, X2
        y     <- 1 + x[,1] - x[,2] + e
        theta <- 1
      } else {
        # Misspecified OM: interaction / nonlinear terms
        y     <- 1 + x[,1] - x[,2] + 0.5 * x[,1] * x[,2] + 0.3 * (x[,2]^2 - 1) + e
        # true expectation can be worked out under Xj ~ N(2,1):
        # E(Y) = 1 + 2 - 2 + 0.5*E[X1*X2] + 0.3*(E[X2^2]-1)
        #      = 1 + 0 + 0.5*4 + 0.3*(5-1) = 1 + 2 + 1.2 = 4.2
        theta <- 4.2
      }
      
      # # Sim new
      # n    <- 1000
      # pcol <- 3
      # x    <- matrix(rnorm(n * pcol, mean = 2, sd = 1), ncol = pcol)
      # # colnames(x) <- c("X1","X2","X3")
      # z    <- x    # no distinction
      # 
      # ## Error term
      # e <- rnorm(n, 0, 1)
      # vx <- rep(1, n)   # keep original definition
      # 
      # ## Outcome model (OR)
      # if (OM) {
      #   # Correct OR: linear in X1, X2
      #   y     <- 1 + x[,1] - x[,2] + e
      #   theta <- 1   # E(Y) = 1 under Xj ~ N(2,1)
      # } else {
      #   # Misspecified OR: interaction / nonlinear terms
      #   y     <- 1 + x[,1] - x[,2] +
      #     0.5 * x[,1] * x[,2] + 0.3 * (x[,2]^2 - 1) + e
      #   # True expectation under Xj ~ N(2,1):
      #   # E(Y) = 1 + 2 - 2 + 0.5*E[X1*X2] + 0.3*(E[X2^2]-1)
      #   #      = 4.2
      #   theta <- 4.2
      # }
      # 
      ## Stratified sampling with fixed sample sizes
      # Define 4 strata based on X2 and X3 (split at 2)
      S2 <- ifelse(x[,2] > 2, 1, 0)
      S3 <- ifelse(x[,3] > 2, 1, 0)
      stratum <- 1 + S2*2 + S3   # (0,0)->1, (0,1)->2, (1,0)->3, (1,1)->4

      # Population stratum sizes
      N_h <- as.integer(tabulate(stratum, nbins = 4))

      # Fixed allocation (equal allocation here; can adjust)
      n_h <- c(150, 100, 100, 50)   # example: oversample stratum 1

      # Draw without replacement from each stratum
      delta <- rep(0, n)
      for (h in 1:4) {
        idx_h <- which(stratum == h)
        samp_h <- sample(idx_h, n_h[h])
        delta[samp_h] <- 1
      }
      
      # Define inclusion probabilities pi
      pi <- rep(NA, n)
      for (h in 1:4) {
        pi[stratum == h] <- n_h[h] / N_h[h]
      }
      
      # ## Response model (RM)
      # if (RM) {
      #   # Correct RM: simple logistic in X1, X2
      #   c  <- -0.5 - 0.25 * x[,2] + 0.5 * x[,3]
      #   pi <- 1 / (1 + exp(-c))
      #   delta = rbinom(n, 1, pi)
      # } else {
      #   # Draw without replacement from each stratum
      #   delta <- rep(0, n)
      #   for (h in 1:4) {
      #     idx_h <- which(stratum == h)
      #     samp_h <- sample(idx_h, n_h[h])
      #     delta[samp_h] <- 1
      #   }
      # 
      #   # Define inclusion probabilities pi
      #   pi <- rep(NA, n)
      #   for (h in 1:4) {
      #     pi[stratum == h] <- n_h[h] / N_h[h]
      #   }
      # }
      
      
      # # Sim setup 2
      # ## Response model (RM)
      # if (RM) {
      #   # Correct RM: logistic in X2, X3 with stronger effects
      #   c  <- -0.7 - 0.6 * x[,2] + 0.4 * x[,3]
      #   pi <- 1 / (1 + exp(-c))
      # } else {
      #   # Misspecified RM: quadratic terms, bigger curvature
      #   c  <- -0.7 + 0.35 * (x[,2] - 2)^2 - 0.25 * (x[,3] - 1.5)^2
      #   pi <- 1 / (1 + exp(-c))
      # }
      # 
      # ## Error term
      # e <- rnorm(n, 0, 1)
      # 
      # ## Outcome model (OM)
      # if (OM) {
      #   y <- 1.5 + 1.2*x[,1] - 0.7*x[,2] + e
      #   theta <- 2.5
      # } else {
      #   y <- 1.5 + 1.2*x[,1] - 0.7*x[,2] + 0.6*x[,1]*x[,2] + 0.4*(x[,2]^2 - 1) + e
      #   theta <- 6.5
      # }
      
      alpha_res = NULL
      theta_res = NULL
      var_res = NULL
      CR_res = NULL
      
      # delta = rbinom(n, 1, pi)
      Index_S = (delta == 1)
      
      y_S = y[Index_S]
      x_S = x[Index_S,]
      z_S = z[Index_S,]
      
      del = NA
      type = "EL"
      
      data = data.frame(y, delta, vx = vx, x = x, z = z)
      data_S = data[Index_S,]
      
      # OM = T; RM = T
      col_RM = 2:3
      # col_RM = 1:pcol
      Rmodel = glm(reformulate(paste0("x.", col_RM), response = "delta"), family = binomial,
                   data = data)
      xR = cbind(1, x[,col_RM])   
      
      # if(RM == T){
      #   Rmodel = glm(reformulate(paste0("x.", col_RM), response = "delta"), family = binomial,
      #                data = data)
      #   # Rmodel = glm(reformulate(c("0", paste0("x.", col_RM)), response = "delta"), family = binomial,
      #   #              data = data)
      #   xR = cbind(1, x[,col_RM])
      # }else if(RM == F){
      #   Rmodel = glm(reformulate(paste0("z.", col_RM), response = "delta"), family = binomial,
      #                data = data)
      #   # Rmodel = glm(delta ~ x.1 + x.2 + z.3 + z.4, family = binomial,
      #   #              data = data) # To be removed
      #   # Rmodel = glm(reformulate(c("0", paste0("z.", col_RM)), response = "delta"), family = binomial,
      #   #              data = data)
      #   xR = cbind(1, z[,col_RM])
      # }
      pihat = predict.glm(Rmodel, data, type = "response")
      # pihat = pi # To be removed
      # pihat = ifelse(pihat > quantile(pihat, 0.05), pihat, quantile(pihat, 0.05))
      
      col_OM = 1:2
      # col_OM = 1:pcol
      Omodel = lm(reformulate(paste0("x.", col_OM), response = "y"), weights = 1 / vx, data = data_S)
      xO = x[, col_OM]; xO_S = xO[Index_S, ]      
      
      # if(OM == T){
      #   Omodel = lm(reformulate(paste0("x.", col_OM), response = "y"), weights = 1 / vx, data = data_S)
      #   xO = x[, col_OM]; xO_S = xO[Index_S, ]
      # }else if(OM == F){
      #   Omodel = lm(reformulate(paste0("z.", col_OM), response = "y"), weights = 1 / vx, data = data_S)
      #   xO = z[, col_OM]; xO_S = xO[Index_S, ]
      # }
      yhat = predict.lm(Omodel, data, type = "response")
      
      safe_rcond <- function(M) {
        if (is.null(M)) return(0)
        if (!all(is.finite(M))) return(0)
        tryCatch(rcond(M), error = function(e) 0)
      }
      
      inv_safe <- function(A, ridge = 1e-8) {
        p <- ncol(A)
        solve(A + diag(ridge, p))
      }
      
      # MLE
      findphi0 = function(phi, xR, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-xR %*% phi)))
        w_phi = delta - pi_phi
        if(returnw){
          return(drop(1 / pi_phi))
        }else{
          return(drop(w_phi %*% Z))
        }
      }
      jacphi0 = function(phi, xR, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-xR %*% phi)))
        return(-t(Z) %*% (xR * (pi_phi * (1 - pi_phi)) ))
      }
      Rmodel0 = glm(delta ~ ., family = binomial,
                    data = data.frame(xR = xR[,-1], delta))
      (nleqslv_res = nleqslv(Rmodel0$coefficients, findphi0, jac = jacphi0, xR = xR, Z = xR, delta = delta,
                             control = list(maxit = 1e5, allowSingular = T), xscalm = "auto",
                             method = "Newton"))
      
      if(nleqslv_res$termcd != 1 & max(abs(findphi0(nleqslv_res$x, xR = xR, Z = xR, delta = delta))) > 1e-5){
        w_S = NA
      }else{
        w_S = findphi0(nleqslv_res$x, xR = xR, Z = xR, delta = delta, returnw = T)
      }
      w_S0 = w_S # MLE
      
      # CBPS
      findphi = function(phi, xR, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-xR %*% phi)))
        w_phi = ifelse(delta, 1 / pi_phi, - 1 / (1 - pi_phi))
        if(returnw){
          return(drop(1 / pi_phi))
        }else{
          return(drop(w_phi %*% Z))
        }
      }
      jacphi = function(phi, xR, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-xR %*% phi)))
        logit_phi = pi_phi / (1 - pi_phi)
        return(-t(Z) %*% (xR * ifelse(delta, 1 / logit_phi, logit_phi) ))
        # return(-t(xR) %*% (Z * ifelse(delta, 1 / logit_phi, logit_phi) ))
      }
      
      Rmodel0 = glm(delta ~ ., family = binomial,
                    data = data.frame(xR = xR[,-1], delta))
      (nleqslv_res = nleqslv(Rmodel0$coefficients, findphi, jac = jacphi, xR = xR, Z = xR, delta = delta,
                             control = list(maxit = 1e5, allowSingular = T), xscalm = "auto",
                             method = "Newton"))
      
      if(nleqslv_res$termcd != 1 & max(abs(findphi(nleqslv_res$x, xR = xR, Z = xR, delta = delta))) > 1e-5){
        w_S = NA
      }else{
        w_S = findphi(nleqslv_res$x, xR = xR, Z = xR, delta = delta, returnw = T)
      }
      w_S1 = w_S # CBPS
      
      # theta_res = c(theta_res, CIPW = sum(y_S * w_S) / sum(w_S)) # CBPS_IPW
      # # theta_res = c(theta_res, CBPS1 = sum(y_S * w_S) / n) # CBPS_IPW
      # theta_res = c(theta_res, CAIPW = (sum(yhat) + sum((y_S - yhat[Index_S]) * w_S)) / n) # CBPS_AIPW
      
      # w = CVXR::Variable(length(d_S))
      # constraints <- list(Z_St %*% w / n == Zbar)
      # 
      # Phi_R <- Maximize(sum(log(w)))
      # 
      # prob <- CVXR::Problem(Phi_R, constraints)
      # # res <- CVXR::solve(prob, solver = "ECOS_BB")
      # res <- CVXR::solve(prob, solver = "SCS")
      # 
      # if(res$status == "optimal" | res$status == "optimal_inaccurate"){
      #   w_S = drop(res$getValue(w))
      # }else{
      #   w_S = NA
      # }
      
      # Tan
      findphi2 = function(phi, xR, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-xR %*% phi)))
        w_phi = ifelse(delta, (1 - pi_phi) / pi_phi, -1)
        if(returnw){
          return(drop(1 / pi_phi))
        }else{
          return(drop(w_phi %*% Z))
        }
      }
      jacphi2 = function(phi, xR, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-xR %*% phi)))
        logit_phi = pi_phi / (1 - pi_phi)
        return(-t(Z) %*% (xR * ifelse(delta, 1 / logit_phi, 0) ))
        # return(-t(xR) %*% (Z * ifelse(delta, 1 / logit_phi, logit_phi) ))
      }
      
      (nleqslv_res = nleqslv(Rmodel0$coefficients, findphi2, jac = jacphi2, xR = xR, Z = xR, delta = delta,
                             control = list(maxit = 1e5, allowSingular = T), xscalm = "auto",
                             method = "Newton"))
      
      if(nleqslv_res$termcd != 1 & max(abs(findphi2(nleqslv_res$x, xR = xR, Z = xR, delta = delta))) > 1e-5){
        w_S = NA
      }else{
        w_S = findphi2(nleqslv_res$x, xR = xR, Z = xR, delta = delta, returnw = T)
      }
      # drop(t(xR[Index_S,]) %*% w_S2[Index_S]); colSums(xR)
      w_S2 = w_S # Tan
      
      # theta_res = c(theta_res, Tan = sum(y_S * w_S) / sum(w_S)) # Tan's regularized calibration (2020)
      
      for(pimethod in c(1)){
        if(pimethod == 1){
          hhat = pihat * xR # MLE
        }else if(pimethod == 2){
          pihat = 1 / w_S1
          hhat = 1 / (1 - pihat) * xR
        }else if(pimethod == 3){
          pihat = 1 / w_S2
          # hhat = pihat / (1 - pihat) * xR # Tan
          hhat = xR # Tan
        } 
        # y_IPW = sum(y_S / pi[Index_S]) / sum(1 / pi[Index_S])
        y_IPW = sum(y_S / pihat[Index_S]) / sum(1 / pihat[Index_S])
        theta_res = c(theta_res, IPW = y_IPW)  # IPW
        if (rcond(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S])) < 1e-12) {
          kappa = rep(NA, ncol(hhat))
        } else {
          kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
                        t(hhat[Index_S,]) %*% ((y_S - y_IPW) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        }
        
        eta = drop(hhat %*% kappa)
        eta[Index_S] = eta[Index_S] + (y_S - y_IPW - eta[Index_S]) / pihat[Index_S]
        var_res = c(var_res, IPW = var(eta) / n)
        # theta_res = c(theta_res, GLS = sum(yhat) / n) # Prediction estimator
        # theta_res = c(theta_res, DR0 = (sum(yhat) / n + sum((y_S - yhat[Index_S]) / pi[Index_S]) / sum(1 / pi[Index_S]))) # AIPW
        # theta_res = c(theta_res, DR = (sum(yhat) + sum((y_S - yhat[Index_S]) / pihat[Index_S])) / sum(1 / pihat[Index_S])) # AIPW
        
        
        
        # theta_res = c(theta_res, DR0 = (sum(yhat) / n + sum((y_S - yhat[Index_S]) / pi[Index_S]) / n)) # AIPW
        theta_res = c(theta_res, AIPW = sum(yhat) / n + sum((y_S - yhat[Index_S]) / pihat[Index_S]) / sum(1 / pihat[Index_S])) # AIPW
        if (rcond(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S])) < 1e-12) {
          gammahat = rep(NA, ncol(hhat[Index_S,]))
        } else {
          kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
                        t(hhat[Index_S,]) %*% ((y_S - yhat[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        }
        
        eta = yhat + drop(hhat %*% kappa)
        eta[Index_S] = eta[Index_S] + (y_S - eta[Index_S]) / pihat[Index_S]
        var_res = c(var_res, AIPW = var(eta) / n)
        
        L_alpha = function(alpha, pihat, n, Index_S, xO, xO_S, y ,yhat, delta){
          nhat = sum(1 / pihat[Index_S])
          q_S = pihat[Index_S]^(alpha - 1)
          Del = colSums(xO_S / pihat[Index_S]) / nhat - colSums(xO) / n
          tmpvec <- tryCatch(drop((xO_S * q_S) %*% solve(t(xO_S * q_S) %*% xO_S / n, Del)) / n * nhat, error = function(e) Inf)
          if(VM){
            return(sum((1 / pihat[Index_S] - tmpvec)^2 * (y[Index_S] - yhat[Index_S])^2))             }else{
            return(sum((1 / pihat[Index_S] - tmpvec)^2))              
          }
        }
        
        # L_alpha <- function(alpha, pihat, n, Index_S, xO, xO_S, y, yhat, delta = delta) {
        #   
        #   # number of covariates
        #   p <- ncol(xO)
        #   
        #   # n̂ = sum δ / pihat
        #   n_hat <- sum(delta / pihat)
        #   
        #   # Δ̂ = weighted mean - population mean
        #   Delta_hat <- (1 / n_hat) * colSums((delta / pihat) * xO) - colMeans(xO)
        #   
        #   # M̂_q(α) = n^{-1} Σ δ q1_i x_i x_i^T, q1 = pihat^{α-1}
        #   q1 <- pihat^(alpha - 1)
        #   Mhat <- matrix(0, p, p)
        #   for (i in seq_len(n)) {
        #     if (delta[i] == 1) {
        #       Mhat <- Mhat + q1[i] * tcrossprod(xO[i, ])
        #     }
        #   }
        #   Mhat <- Mhat / n
        #   Minv <- inv_safe(Mhat)
        #   
        #   # L̂(α) = n^{-1} Σ δ { n/n̂ * pihat_i - Δ̂^T M^{-1} q1_i x_i }^2
        #   tvals <- numeric(n)
        #   for (i in seq_len(n)) {
        #     if (delta[i] == 1) {
        #       tvals[i] <- (n / n_hat) * pihat[i] -
        #         as.numeric(crossprod(Delta_hat, Minv %*% (q1[i] * xO[i, ])))
        #     } else {
        #       tvals[i] <- 0
        #     }
        #   }
        #   
        #   mean(delta * tvals^2)
        # }
        
        
        alphahat = optimize(L_alpha, pihat = pihat, n = n, Index_S = Index_S, xO = cbind(1, xO), xO_S= cbind(1, xO_S), y = y, yhat = yhat,
                            delta = delta,
                            interval = c(-1, 3))$minimum        
        
        q = pihat^(alphahat - 1); q_S = q[Index_S]
        
        yhatq = drop(cbind(1, xO) %*% solve(t(cbind(1, xO_S) * q_S) %*% cbind(1, xO_S), t(cbind(1, xO_S) * q_S) %*% y_S))
        
        # theta_res = c(theta_res, DR0 = (sum(yhat) / n + sum((y_S - yhat[Index_S]) / pi[Index_S]) / n)) # AIPW
        theta_res = c(theta_res, AIPW2 = sum(yhatq) / n + sum((y_S - yhatq[Index_S]) / pihat[Index_S]) / sum(1 / pihat[Index_S])) # AIPW
        if (rcond(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S])) < 1e-12) {
          gammahat = rep(NA, ncol(hhat[Index_S,]))
        } else {
          kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
                        t(hhat[Index_S,]) %*% ((y_S - yhatq[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        }
        
        eta = yhatq + drop(hhat %*% kappa)
        eta[Index_S] = eta[Index_S] + (y_S - eta[Index_S]) / pihat[Index_S]
        var_res = c(var_res, AIPW2 = var(eta) / n)
        
        
        gprime1 = function(x, type, del){ # Inverse of gprime(d)
          switch(type,
                 SL = rep(1, length(x)),
                 EL = x^2,
                 ET = x,
                 CE = x * (x - 1),
                 HD = 2 * x^(1.5),
                 Hb = ifelse(abs(x) < del, 1, NA), # ???
                 PH = (1 + (x / del)^2)^(1.5),
                 GE = x^(1 - del))
        }
        
        G = function(x, type, del){
          switch(type,
                 SL = x^2/2,
                 EL = -log(x),
                 ET = x * (log(x) - 1),
                 CE = (x-1) * log(x-1) - x * log(x),
                 GE = x^(del + 1) / del / (del + 1),
                 HD = -4 * sqrt(x),
                 PH = del^2 * sqrt(1 + (x / del)^2))
        }
        
        g <- function(pihat, type, del = NULL) {
          if (type == "EL") {
            u_vec_tmp <- -pihat
          } else if (type == "ET") {
            u_vec_tmp <- -log(pihat)
          } else if (type == "GE") {
            u_vec_tmp <- (1 / pihat)^del / del
          } else if (type == "SL") {
            u_vec_tmp <- 1 / pihat
          } else {
            stop("Unknown type. Must be one of: 'EL', 'ET', 'GE', 'SL'.")
          }
          return(u_vec_tmp)
        }
        
        f = function(lambda, d_S, v_S, Z_S, Zbar, type, del, n, ..., returnw = F){
          # print(lambda)
          # w_S = d_S * ginv(drop(Z_S %*% lambda), type = type)
          if(type == "HD" & any(Z_S %*% lambda / v_S >= 0)) return(rep(Inf, length(lambda)))
          if(!returnw & type == "GE" & any(drop(Z_S %*% lambda / v_S) * del <= 0)) return(rep(Inf, length(lambda)))
          if(type == "PH" & any(abs(Z_S %*% lambda / v_S) >= del)) return(rep(Inf, length(lambda)))
          
          if(type == "SL"){
            w_S = d_S * drop(Z_S %*% lambda / v_S)
          }else if(type == "EL"){
            w_S = -d_S / drop(Z_S %*% lambda / v_S)
          }else if(type == "ET"){
            w_S = d_S * exp(drop(Z_S %*% lambda / v_S))
          }else if(type == "CE"){
            w_S = d_S / (1 - exp(drop(Z_S %*% lambda / v_S)))
          }else if(type == "HD"){
            w_S = d_S / drop(Z_S %*% lambda / v_S)^2
          }else if(type == "PH"){
            w_S = d_S / sqrt(1 / drop(Z_S %*% lambda / v_S)^2 - 1 / del^2)
          }else if(type == "GE"){
            w_S = d_S * (drop(Z_S %*% lambda / v_S) * del)^(1 / del)
          }
          # print(w_S)
          if(!returnw & type != "SL" & any(w_S <= 0)) return(rep(Inf, length(lambda)))
          if(type == "CE" & any(w_S <= 1)) return(rep(Inf, length(lambda)))
          if(returnw == T){
            return(w_S)
          }else{
            return(colSums(Z_S * w_S) / n - Zbar)
          }
        }
        
        h = function(lambda, d_S, v_S, Z_S, Z_St, Zbar, type, del, n){
          # return(Z_St %*% (Z_S * d_S * ginvprime(drop(Z_S %*% lambda), type = type) / n))
          if(type == "SL"){
            return(Z_St %*% (Z_S / v_S * d_S) / n)
          }else if(type == "EL"){
            w_S = -d_S / drop(Z_S %*% lambda / v_S)
            return(Z_St %*% (Z_S / v_S * (w_S^2 / d_S)) / n)
          }else if(type == "ET"){
            w_S = d_S * exp(drop(Z_S %*% lambda / v_S))
            return(Z_St %*% (Z_S / v_S * w_S) / n)
          }else if(type == "CE"){
            p_Stmp = 1 / (1 - exp(drop(Z_S %*% lambda / v_S)))
            return(Z_St %*% (Z_S / v_S * d_S * p_Stmp * (p_Stmp - 1)) / n)
          }else if(type == "HD"){
            return(Z_St %*% (Z_S / v_S * (-2 * d_S / drop(Z_S %*% lambda / v_S)^3)) / n)
          }else if(type == "PH"){
            return(Z_St %*% (Z_S / v_S * (d_S * (1 - (drop(Z_S %*% lambda / v_S) / del)^2)^(-1.5))) / n)
          }else if(type == "GE"){
            return(Z_St %*% (Z_S / v_S * (d_S * drop(Z_S %*% lambda  / v_S * del)^(1 / del - 1)) / n))
          }
        }
        
        targetftn0 = function(W, d_S, Z_S, Z_St, init, Zbar, type, del,..., returnw = F){
          if(type == "ET" | type == "SL" | type == "PH"){
            if(W < 0) return(.Machine$double.xmax)
          }else if(type == "EL" | type == "CE" | type == "HD"){
            if(W > 0) return(.Machine$double.xmax)
          }
          alphaHT = Zbar[length(Zbar)]
          Zbar[length(Zbar)] <- W
          nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, v_S = v_S,
                                Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n,
                                method = "Newton", control = list(maxit = 1e5, allowSingular = T))
          # if(nleqslv_res$termcd != 1 & nleqslv_res$termcd != 2){
          if(nleqslv_res$termcd != 1){
            if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n))) > 1e-5)
              return(.Machine$double.xmax)
          }
          w_S = f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n, returnw = T)
          # if(any(is.infinite(w_S))) return(.Machine$double.xmax)
          
          
          
          if(returnw == F){
            return(sum(G(w_S, type = type, del = del)) - n * W)
            # return(sum(G(w_S, type = type, del = del)) - n * (alphaHT + 1) * log(abs(W + 1)))
          }else{
            return(w_S)
          }
        }
        
        # d_S = rep(1, sum(delta)); u_vec = 1 / pihat # Reg
        # # d_S = vx[Index_S]; u_vec = -pihat * vx # EL2
        # u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec); 
        # Z_S = cbind(1, xO_S, u_vec_S); Zbar = c(1, colMeans(xO), Uhat); Z_St = t(Z_S)
        # 
        # init = rep(0, length(Zbar)); init[length(init)] = 1
        # 
        # nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
        #                       Z_St = Z_St, Zbar = Zbar, type = "SL", del = del,
        #                       method = "Newton", control = list(maxit = 1e5, allowSingular = T),
        #                       xscalm = "auto")
        # if(nleqslv_res$termcd != 1){
        #   if(max(abs(f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = "SL", del = del))) > 1e-5)
        #     w_S = NA
        # }else{
        #   w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = "SL", del = del, returnw = T)             
        # }
        # 
        # theta_res = c(theta_res, Reg = sum(y_S * w_S) / n) # Reg
        
        # type = "GE"; del = -0.5
        # type = "EL"
        # type = "ET"

        DG = function(w_S0, w_S, q_S, type, del){
          sum((G(w_S, type, del) - G(w_S0, type, del) -
          (g(1 / w_S0, type, del) * (w_S - w_S0))) / q_S)
        }
        
        for(type in c("EL", "ET", "GE", "SL")){
          if(type == "GE") del = -0.5
          d_S = rep(1, sum(delta)); 
          vx = rep(1, length(delta)); v_S = vx[Index_S];          
          u_vec_tmp = g(pihat, type, del)
          
          L_alpha <- function(alpha, pihat, n, Index_S, xO, xO_S, y, yhat, u_vec_tmp, d_S, type, del){
            q  <- pihat^(alpha - 1)
            vx <- 1 / q; v_S <- vx[Index_S]
            u_vec <- matrix(u_vec_tmp * vx, ncol = 1)
            u_vec_S <- u_vec[Index_S, , drop = FALSE]; Uhat <- colMeans(u_vec)
            Z_S  <- cbind(1, xO_S, u_vec_S); Zbar <- c(1, colMeans(xO), Uhat); Z_St <- t(Z_S)
            init <- rep(0, length(Zbar)); init[length(init)] <- 1
            
            nleqslv_res <- nleqslv(init, f, jac = h,
                                   d_S = d_S, v_S = v_S, Z_S = Z_S, Z_St = Z_St,
                                   Zbar = Zbar, type = type, del = del, n = n,
                                   method = "Newton",
                                   control = list(maxit = 1e5, allowSingular = TRUE),
                                   xscalm = "auto")
            if (nleqslv_res$termcd != 1) {
              if (max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S,
                            Zbar = Zbar, type = type, del = del, n = n))) > 1e-5) {
                w_S <- NA
              }
            } else {
              w_S <- f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar,
                       type = type, del = del, n = n, returnw = TRUE)
            }
            
            Gvals <- gprime1(w_S, type = type, del = del)
            Amat  <- Z_St %*% (Z_S * Gvals / d_S / v_S)
            if (safe_rcond(Amat) < 1e-12) {
              gammahat <- rep(NA_real_, ncol(Z_S))
            } else {
              gammahat <- solve(Amat, Z_St %*% (y_S * Gvals / d_S / v_S))
            }
            yhat2 <- drop(cbind(1, xO, u_vec) %*% gammahat)
            
            if (VM) {
              return(sum(w_S^2 * (y[Index_S] - yhat2[Index_S])^2))
            } else {
              return(sum(w_S^2))
            }
          }
          
          ## 1) get alphahat2 and q once (used in cases 2 & 4)
          alphahat2 <- optimize(
            L_alpha,
            interval = c(-1, 3),
            pihat = pihat, n = n, Index_S = Index_S, xO = xO, xO_S = xO_S,
            y = y, yhat = NULL, u_vec_tmp = u_vec_tmp, d_S = d_S, type = type, del = del
          )$minimum
          q <- pihat^(alphahat2 - 1)
          
          ## 2) define the 4 cases’ controls
          case_id  <- 1:4
          vx_list  <- list(
            rep(1, length(delta)),    # case 1: constant weighting
            1 / q,                    # case 2: q-weighting
            rep(1, length(delta)),    # case 3: constant + orthogonal
            1 / q                     # case 4: q-weighting + orthogonal
          )
          ortho_list <- c(FALSE, FALSE, TRUE, TRUE)
          DGq_list   <- list(
            1,                        # case 1
            q[Index_S],               # case 2
            1,                        # case 3
            q[Index_S]                # case 4
          )
          
          ## 3) loop over cases
          for (ii in seq_along(case_id)) {
            cid <- case_id[[ii]]
            vx  <- vx_list[[ii]]
            do_ortho <- ortho_list[[ii]]
            DG_q <- DGq_list[[ii]]
            
            v_S   <- vx[Index_S]
            u_vec <- matrix(u_vec_tmp * vx, ncol = 1)
            
            if (do_ortho) {
              u_vec <- cbind(
                -xR * ((1 - pihat) / pihat / gprime1(1 / pihat, type = type, del = del) * vx),
                u_vec
              )
            }
            
            u_vec_S <- u_vec[Index_S, , drop = FALSE]
            Uhat    <- colMeans(u_vec)
            Z_S     <- cbind(1, xO_S, u_vec_S)
            Zbar    <- c(1, colMeans(xO), Uhat)
            Z_St    <- t(Z_S)
            Z       <- cbind(1, xO, u_vec)
            
            init <- rep(0, length(Zbar)); init[length(init)] <- 1
            
            nleqslv_res <- nleqslv(init, f, jac = h,
                                   d_S = d_S, v_S = v_S, Z_S = Z_S, Z_St = Z_St,
                                   Zbar = Zbar, type = type, del = del, n = n,
                                   method = "Newton",
                                   control = list(maxit = 1e5, allowSingular = TRUE),
                                   xscalm = "auto")
            
            if (nleqslv_res$termcd != 1) {
              if (max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S,
                            Zbar = Zbar, type = type, del = del, n = n))) > 1e-5) {
                w_S <- NA
              }
            } else {
              w_S <- f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar,
                       type = type, del = del, n = n, returnw = TRUE)
            }
            
            ## theta
            theta_res <- c(theta_res, setNames(sum(y_S * w_S) / n, paste(type, cid, sep = "")))
            
            ## gammahat and yhat2
            Gvals <- gprime1(w_S, type = type, del = del)
            Amat  <- Z_St %*% (Z_S * Gvals / d_S / v_S)
            if (safe_rcond(Amat) < 1e-12) {
              gammahat <- rep((NA_real_), ncol(Z_S))
            } else {
              gammahat <- solve(Amat, Z_St %*% (y_S * Gvals / d_S / v_S))
            }
            yhat2 <- drop(Z %*% gammahat)
            
            ## kappa (recompute per case)
            HH <- t(hhat[Index_S, ]) %*% (hhat[Index_S, ] * (1 / pihat[Index_S] - 1) / pihat[Index_S])
            if (rcond(HH) < 1e-12) {
              kappa <- rep(NA_real_, ncol(hhat))
            } else {
              kappa <- solve(HH,
                             t(hhat[Index_S, ]) %*%
                               ((y_S - yhat2[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
            }
            
            ## influence/variance
            eta <- yhat2 + drop(hhat %*% kappa)
            if (do_ortho) {
              # orthogonal calibration cases (3 & 4)
              eta[Index_S] <- eta[Index_S] + (y_S - yhat2[Index_S]) * w_S
            } else {
              # non-orthogonal cases (1 & 2)
              eta[Index_S] <- eta[Index_S] + (y_S - yhat2[Index_S]) * w_S -
                drop(hhat %*% kappa)[Index_S] / pihat[Index_S]
            }
            var_res <- c(var_res, setNames(var(eta) / n, paste(type, cid, sep = "")))
            
            # ## alpha metric
            # alpha_res <- c(alpha_res, setNames(
            #   DG(1 / pihat[Index_S], w_S, DG_q, type, del),
            #   paste(type, cid, sep = "")
            # ))
            
            if(ii == 1){
              w_S0 <- w_S
            }else if(ii == 2){
              w_S1 <- w_S
            }
            if(ii == 3){
              alpha_res <- c(alpha_res, setNames(
                DG(w_S0, w_S, DG_q, type, del),
                paste(type, cid, sep = "")
              ))
            }else if(ii == 4){
              alpha_res <- c(alpha_res, setNames(
                DG(w_S1, w_S, DG_q, type, del),
                paste(type, cid, sep = "")
              ))
            }else{
              alpha_res <- c(alpha_res, setNames(
                DG(1 / pihat[Index_S], w_S, DG_q, type, del),
                paste(type, cid, sep = "")
              ))
            }
          }
          
        }
        
        
        CR_res = ifelse(abs(theta_res - theta) < qnorm(0.975) * sqrt(var_res), 1, 0)
      }
      # nlmres = NULL
      # alpha_vec = c(alpha_vec, alphahat2)
      alpha_vec = c(alpha_vec, alpha_res)
      theta_mat = cbind(theta_mat, theta_res - theta)
      var_mat = cbind(var_mat, var_res)
      CR_mat = cbind(CR_mat, CR_res)
    }
    list(theta_mat, alpha = alpha_vec, var_mat, CR_mat)
  }
# final_res
# save.image(paste("rdata/", timenow0, ".RData", sep = ""))
final_res1 = lapply(final_res, function(x) x[[1]])

stopCluster(cl)
timenow2 = Sys.time()
print(timenow2 - timenow1)
# if(sum(!sapply(final_res1, function(x) is.numeric(unlist(x)))) != 0) stop(paste(final_res1))
paste("# of failure:", sum(!sapply(final_res1, function(x) is.numeric(unlist(x))))); final_res0 = final_res1
print(final_res0[!sapply(final_res0, function(x) is.numeric(unlist(x)))])
final_res1 = final_res1[sapply(final_res1, function(x) is.numeric(unlist(x)))]

# res1 <- Reduce(`+`, final_res1) / length(final_res1)
# colnames(res1) = rnames

# BIAS
res1 <- matrix(rowMeans(do.call(cbind, lapply(final_res1, c)), na.rm = TRUE), 
               ncol = ncol(final_res1[[1]]))
colnames(res1) = rnames; row.names(res1) = row.names(final_res1[[1]])

# RMSE
final_res2 = lapply(final_res1, function(x) x^2)
final_res2 = final_res2[sapply(final_res1, function(x) is.numeric(unlist(x)))]
# res2 <- sqrt(Reduce(`+`, final_res2) / length(final_res2))
# colnames(res2) = rnames

# boxplot(t(sapply(final_res1, function(x) x[c(-1,-3),3])))
# summary(sapply(final_res1, function(x) x[6,3]))

res2 <- sqrt(matrix(rowMeans(do.call(cbind, lapply(final_res2, c)), na.rm = TRUE), 
                    ncol = ncol(final_res2[[1]])))
colnames(res2) = rnames; row.names(res2) = row.names(final_res1[[1]])

# colSums(apply(do.call(cbind, lapply(final_res1, c)), 1, is.na))
na_res = matrix(colSums(apply(do.call(cbind, lapply(final_res2, c)), 1, is.na)),
                ncol = ncol(final_res2[[1]]))
colnames(na_res) = rnames; row.names(na_res) = row.names(final_res1[[1]])
na_res

res3 <- matrix(apply(do.call(cbind, lapply(final_res1, c)), 1, 
                     function(x) sqrt(var(x, na.rm = TRUE) * (length(x)-1) / (length(x)))), 
               ncol = ncol(final_res1[[1]]))
colnames(res3) = rnames; row.names(res3) = row.names(final_res1[[1]])

round(res1, 4)
round(res2, 4)
# round(res2^2, 6)

# xtable::xtable(res1, digits = 3)
# xtable::xtable(res2, digits = 3)

Ncol = ncol(final_res2[[1]])
Nrow = nrow(final_res2[[1]])

for(ncols in 1:Ncol){
  boxres = t(do.call(cbind, lapply(final_res1, c)))[, ((ncols - 1) * Nrow + 1) : (ncols * Nrow)]
  colnames(boxres) <- row.names(final_res1[[1]])
  if(!interactive()) png(paste("boxplot", ncols, ".png", sep = ""))
  boxplot(boxres, main = paste(paste(c("RP", "OR"), ifelse(modelcases[ncols,], "C", "M"), sep = ":"), collapse = ", "))
  abline(h = 0, col = "red")
  if(!interactive()) dev.off()
}


# if(!interactive()) png("boxplots.png", width = 960, height = 960)
# # op <- par(); par(mfrow = c(2,2), mar = rep(3, 4))
# for(ncols in 1:Ncol){
#   boxres = t(do.call(cbind, lapply(final_res1, c)))[, ((ncols - 1) * Nrow + 1) : (ncols * Nrow)] + theta
#   colnames(boxres) <- row.names(final_res1[[1]])
#   boxplot(boxres, main = paste(paste(c("RP", "OR"), ifelse(modelcases[ncols,], "C", "M"), sep = ":"), collapse = ", "), 
#           ylim = c(theta - 10, theta + 10), cex.axis=1.2, cex.main = 1.2)
#   # boxplot(boxres, main = rnames[ncols], ylim = c(theta - 10, theta + 10), cex.axis=0.8)
#   abline(h = theta, col = "red")
# }
# # par(op)
# if(!interactive()) dev.off()

final_res_alpha = lapply(final_res, function(x) x[[2]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
alpha_df <- as.data.frame(do.call("rbind", final_res_alpha))
alphanames = names(final_res_alpha[[1]])
colnames(alpha_df) <- alphanames

colnames(alpha_df)[9:12] <- c("HD1", "HD2", "HD3", "HD3")


par(mfrow = c(2,2))

# for(cnt in 1:4){
#   idx = c(1,5,9,13)[cnt]
#   # tmp = c(2, 1/2, 1, 1/8)[cnt] 
#   df = 3
#   qqplot(qchisq(ppoints(length(alpha_df[,idx])), df = df), alpha_df[,idx],
#          main = paste("Q-Q Plot for Chi-squared(", df, ') for', alphanames[idx]),
#          xlab = "Theoretical Quantiles",
#          ylab = "Sample Quantiles")   
#   abline(0,1)
# }

for(cnt in 1:4){
  # cnt = 3
  idx = c(1,5,9,13)[cnt]
  # tmp = c(2, 1/2, 1, 1/8)[cnt] 
  tmpdat = reshape2::melt(alpha_df[order(rowSums(alpha_df[,c(idx,idx+2)])),c(idx,idx+2)])
  tmpdat = cbind(tmpdat, index = c(1:nrow(alpha_df), 1:nrow(alpha_df)))
  tmpdat$variable <- factor(tmpdat$variable, levels = rev(levels(tmpdat$variable)))
  tmptext = paste("Bregman Divergence for", c("EL", "ET", "HD", "SL")[cnt])
  p <- ggplot(tmpdat, aes(x = factor(index), y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    theme_minimal() +
    labs(
      title = tmptext,
      x = "Index",
      y = "Divergence"
    ) + theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
  print(p)
  
  tmpbd = alpha_df[order(rowSums(alpha_df[,c(idx,idx+2)])),c(idx,idx+2)]
  tmpbd = cbind(tmpbd, rowSums(tmpbd))
  print(mean(tmpbd[,1] / tmpbd[,3]))
  
}



head(alpha_df[order(rowSums(alpha_df[,c(idx,idx+2)])),c(idx,idx+2)])

# dev.off()

# alpha_df_long <- pivot_longer(alpha_df, cols = everything(), names_to = "Column", values_to = "Value")
# alpha_df_long$Column <- factor(alpha_df_long$Column, levels = rnames)

# # Plot using ggplot2
# pGG <- ggplot(alpha_df_long, aes(x = Value, fill = Column)) +
#   geom_density(alpha = 0.5) +
#   scale_fill_brewer(palette = "Pastel1") + # Optional: use a color palette
#   theme_minimal() + # Optional: use a minimal theme for the plot
#   labs(title = "Density plot of alpha", x = "alpha", y = "Density") +
#   theme(legend.title = element_blank()) # Optional: remove legend title
# pGG

if(!interactive()) ggsave("hist_alpha.png", plot = pGG)

gsub("\\\\addlinespace", "", kable(cbind(res1, res2) * 10^3, "latex", booktabs = TRUE, digits = 1) %>%
       kable_styling())

# gsub("\\\\addlinespace", "", kable(res1, "latex", booktabs = TRUE, digits = 3) %>%
#   kable_styling())
# 
# gsub("\\\\addlinespace", "", kable(res2, "latex", booktabs = TRUE, digits = 3) %>%
#        kable_styling())

final_res_var = lapply(final_res, function(x) x[[3]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
res_var <- matrix(rowMeans(do.call(cbind, lapply(final_res_var, c)), na.rm = TRUE), 
                  ncol = ncol(final_res_var[[3]]))
colnames(res_var) = rnames; row.names(res_var) = row.names(final_res_var[[3]])

na_res = matrix(colSums(apply(do.call(cbind, lapply(final_res_var, c)), 1, is.na)),
                ncol = ncol(final_res_var[[1]]))
colnames(na_res) = rnames; row.names(na_res) = row.names(final_res1[[1]])
na_res

round((res_var - res3^2) / res3^2, 3)

final_res_CR = lapply(final_res, function(x) x[[4]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
res_CR <- matrix(rowMeans(do.call(cbind, lapply(final_res_CR, c)), na.rm = TRUE), 
                 ncol = ncol(final_res_CR[[4]]))
colnames(res_CR) = rnames; row.names(res_CR) = row.names(final_res_CR[[4]])
round(res_CR, 3)

gsub("\\\\addlinespace", "", kable(cbind((res_var - res3^2) / res3^2, res_CR) * 100, "latex", booktabs = TRUE, digits = 1) %>%
       kable_styling())

# gsub("\\\\addlinespace", "", kable((res_var - res3^2) / res3^2, "latex", booktabs = TRUE, digits = 3) %>%
#        kable_styling())
# 
# gsub("\\\\addlinespace", "", kable(res_CR, "latex", booktabs = TRUE, digits = 3) %>%
#        kable_styling())


