# Simulation setup

if(!interactive()){
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
}else{
  args <- c(50)
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# install.packages( "MatrixModels", type="win.binary" )  

library(nleqslv)
library(CVXR)
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
library(caret)
library(tidyverse)
library(kableExtra)

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

# modelcases = expand.grid(c(T,F),c(T,F),c(T,F))
modelcases = expand.grid(c(T,F),c(T,F))
rnames <- apply(apply(modelcases, 2, ifelse, "C", "M"), 1, paste, collapse = "")

theta = 210

final_res <- foreach(
  simnum = 1:SIMNUM, 
  .packages = c("nleqslv", "CVXR", "caret"),
  .errorhandling="pass") %dopar% {
    theta_mat = NULL
    var_mat = NULL
    CR_mat = NULL
    alpha_vec = NULL
    for(cnt in 1:nrow(modelcases)){
      RM = modelcases[cnt, 1]
      OM = modelcases[cnt, 2]
      # VM = modelcases[cnt, 3]
      
      # Simulation setup of
      # Kang & Schafer, 2007
      n = 1000
      pcol = 4; # pcol = 4 is True
      x = matrix(rnorm(n * pcol, 0, 1), nc= pcol)
      z = cbind(exp(x[,1] /  2), x[,2] / (1 + exp(x[,1])) + 10,
                (x[,1] * x[,3] / 25 + 0.6)^3, (x[,2] + x[,4] + 20)^2) # True
      
      # z = cbind(exp(x[,1] /  2), 1 / (1 + exp(x[,2])) + 10,
      #           (x[,1] * x[,3] / 25 + 0.6)^3, (x[,2] + x[,4] + 20)^2) # To be removed
      
      # z = cbind((x[,1] + 1)^2, (x[,2] + 1)^2, (x[,3] + 1)^2, (x[,4] + 1)^2) # To be removed
      
      # pi = 1 / (1 + exp(-(-0.5 * x[,1] + 0.25 * x[,2] - 3)))
      
      # pi = 1 / (1 + exp(-(-x[,1] + 0.5 * x[,2] -
      #                       0.25 * x[,3] - 0.1 * x[,4] - 2.5)))
      
      pi = 1 / (1 + exp(-(-x[,1] + 0.5 * x[,2] -
                            0.25 * x[,3] - 0.1 * x[,4]))) # True
      
      # pi = ifelse(pt(-rnorm(n, 0, 1) - 2, 3) >.7, .7, pi) # To be removed
      
      # pi = runif(n, 0, 1) # To be removed
      
      # VM = T
      # if(VM == T){
      #   vx = rep(160, n)
      # }else{
      vx = (x[,1]^2 + exp(x[,3]) / 1.65) * 80
      #   # vx = ifelse(x[,1]^2 > 2, 2, x[,1]^2)
      # }
      # vx = rep(500, n)
      # vx = rep(160, n)
      # vx = rep(1, n) # To be removed
      e = rnorm(n, 0, sqrt(vx))
      
      y = theta + 27.4 * x[,1] + 13.7 * x[,2] +
        13.7 * x[,3]+ 13.7 * x[,4] + e
      
      # sqrt((27.4^2 + 13.7^2 * 3 + mean(vx / pi)) / n)
      
      # pi = 1 / (1 + exp(-(-x[,1] + 0.5 * x[,2] -
      #                       0.25 * x[,3] - 0.1 * x[,4] - e / 40))) # To be corrected
      
      # theta = mean(y) # To be corrected
      # y = theta + 27.4 * x[,1] + 13.7 * x[,2] +
      #   13.7 * x[,3]+ 13.7 * x[,4] + e + 50*(sqrt(pi) - 2/ 3) # To be removed
      # sd(y-e) / sd(e)
      # plot(y-e, y)
      
      
      # summary(z)
      # sd(y-e) / sd(e)
      # mean(y)
      # sum(y * pi) / sum(pi) # Treated mean
      # sum(y * (1 - pi)) / sum(1 - pi)
      
      theta_res = NULL
      var_res = NULL
      CR_res = NULL
      
      delta = rbinom(n, 1, pi)
      Index_S = (delta == 1)
      
      y_S = y[Index_S]
      x_S = x[Index_S,]
      z_S = z[Index_S,]
      
      del = NA
      type = "EL"
      
      data = data.frame(y, delta, vx = vx, x = x, z = z)
      data_S = data[Index_S,]
      
      # OM = T; RM = T
      if(RM == T){
        Rmodel = glm(reformulate(paste0("x.", 1:pcol), response = "delta"), family = binomial,
                     data = data)
        # Rmodel = glm(reformulate(c("0", paste0("x.", 1:pcol)), response = "delta"), family = binomial,
        #              data = data)
        x_S0 = cbind(1, x_S); x0 = cbind(1, x)
      }else{
        Rmodel = glm(reformulate(paste0("z.", 1:pcol), response = "delta"), family = binomial,
                     data = data)
        # Rmodel = glm(delta ~ x.1 + x.2 + z.3 + z.4, family = binomial,
        #              data = data) # To be removed
        # Rmodel = glm(reformulate(c("0", paste0("z.", 1:pcol)), response = "delta"), family = binomial,
        #              data = data)
        x_S0 = cbind(1, z_S); x0 = cbind(1, z)
      }
      pihat = predict.glm(Rmodel, data, type = "response")
      # pihat = pi # To be removed
      # pihat = ifelse(pihat > quantile(pihat, 0.05), pihat, quantile(pihat, 0.05))
      
      if(OM == T){
        Omodel = lm(reformulate(paste0("x.", 1:pcol), response = "y"), weights = 1 / vx, data = data_S)
      }else{
        Omodel = lm(reformulate(paste0("z.", 1:pcol), response = "y"), weights = 1 / vx, data = data_S)
      }
      yhat = predict.lm(Omodel, data, type = "response")
      
      # x_S0 = x_S; x0 = x
      # z = scale(z); attributes(z)[2:3] <- NULL; z_S = z[Index_S,];
      if(OM == F){
        x_S = z_S; x = z
        # x_S[,c(3,4)] = z_S[,c(3,4)]; x[,c(3,4)] = z[,c(3,4)] # To be removed
      }
      
      findphi = function(phi, x0, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-x0 %*% phi)))
        w_phi = ifelse(delta, 1 / pi_phi, - 1 / (1 - pi_phi))
        if(returnw){
          return(drop(1 / pi_phi))
        }else{
          return(drop(w_phi %*% Z))
        }
      }
      jacphi = function(phi, x0, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-x0 %*% phi)))
        logit_phi = pi_phi / (1 - pi_phi)
        return(-t(Z) %*% (x0 * ifelse(delta, 1 / logit_phi, logit_phi) ))
        # return(-t(x0) %*% (Z * ifelse(delta, 1 / logit_phi, logit_phi) ))
      }
      
      Rmodel0 = glm(reformulate(paste0("x0.", 1:pcol), response = "delta"), family = binomial,
                    data = data.frame(x0 = x0[,-1], delta))
      
      (nleqslv_res = nleqslv(Rmodel0$coefficients, findphi, jac = jacphi, x0 = x0, Z = x0, delta = delta,
                             control = list(maxit = 1e5, allowSingular = T), xscalm = "auto",
                             method = "Newton"))
      
      if(nleqslv_res$termcd != 1 & max(abs(findphi(nleqslv_res$x, x0 = x0, Z = x0, delta = delta))) > 1e-5){
        w_S = NA
      }else{
        w_S = findphi(nleqslv_res$x, x0 = x0, Z = x0, delta = delta, returnw = T)
      }
      w_S1 = w_S
      
      # theta_res = c(theta_res, CSIPW = sum(y_S * w_S) / sum(w_S)) # CBPS_IPW
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
      
      findphi2 = function(phi, x0, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-x0 %*% phi)))
        w_phi = ifelse(delta, (1 - pi_phi) / pi_phi, -1)
        if(returnw){
          return(drop(1 / pi_phi))
        }else{
          return(drop(w_phi %*% Z))
        }
      }
      jacphi2 = function(phi, x0, Z, delta, ..., returnw = F){
        pi_phi = drop(1 / (1 + exp(-x0 %*% phi)))
        logit_phi = pi_phi / (1 - pi_phi)
        return(-t(Z) %*% (x0 * ifelse(delta, 1 / logit_phi, 0) ))
        # return(-t(x0) %*% (Z * ifelse(delta, 1 / logit_phi, logit_phi) ))
      }
      
      (nleqslv_res = nleqslv(Rmodel0$coefficients, findphi2, jac = jacphi2, x0 = x0, Z = x0, delta = delta,
                             control = list(maxit = 1e5, allowSingular = T), xscalm = "auto",
                             method = "Newton"))
      
      if(nleqslv_res$termcd != 1 & max(abs(findphi2(nleqslv_res$x, x0 = x0, Z = x0, delta = delta))) > 1e-5){
        w_S = NA
      }else{
        w_S = findphi2(nleqslv_res$x, x0 = x0, Z = x0, delta = delta, returnw = T)
      }
      # drop(t(x0[Index_S,]) %*% w_S2[Index_S]); colSums(x0)
      w_S2 = w_S
      
      # theta_res = c(theta_res, Tan = sum(y_S * w_S) / sum(w_S)) # Tan's regularized calibration (2020)
      
      for(pimethod in 1){
        if(pimethod == 1){
          hhat = pihat * x0
        }else if(pimethod == 2){
          pihat = 1 / w_S1
          hhat = 1 / (1 - pihat) * x0
        }else if(pimethod == 3){
          pihat = 1 / w_S2
          hhat = pihat / (1 - pihat) * x0
        } 
        
        theta_res = c(theta_res, SIPW = sum(y_S / pihat[Index_S]) / sum(1 / pihat[Index_S]))  # SIPW
        var_res = c(var_res, SIPW = 0)
        # theta_res = c(theta_res, GLS = sum(yhat) / n) # Prediction estimator
        # theta_res = c(theta_res, DR0 = (sum(yhat) / n + sum((y_S - yhat[Index_S]) / pi[Index_S]) / sum(1 / pi[Index_S]))) # AIPW
        # theta_res = c(theta_res, DR = (sum(yhat) + sum((y_S - yhat[Index_S]) / pihat[Index_S])) / sum(1 / pihat[Index_S])) # AIPW
        
        
        
        # theta_res = c(theta_res, DR0 = (sum(yhat) / n + sum((y_S - yhat[Index_S]) / pi[Index_S]) / n)) # AIPW
        theta_res = c(theta_res, AIPW = (sum(yhat) + sum((y_S - yhat[Index_S]) / pihat[Index_S])) / n) # AIPW
        
        kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
                      t(hhat[Index_S,]) %*% ((y_S - yhat[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        eta = yhat + drop(hhat %*% kappa)
        eta[Index_S] = eta[Index_S] + (y_S - yhat[Index_S] - drop(hhat %*% kappa)[Index_S]) / pihat[Index_S]
        var_res = c(var_res, AIPW = var(eta) / n)
        
        # Z_S = cbind(1, x_S); Zbar = c(1, colMeans(x));  Z_St = t(Z_S); d_S = rep(1, sum(delta)); 
        
        
        
        
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
                 HD = -2 * sqrt(x),
                 PH = del^2 * sqrt(1 + (x / del)^2))
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
        # Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat); Z_St = t(Z_S)
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
        
        d_S = rep(1, sum(delta)); 
        v_S = rep(1, sum(delta)); u_vec = -pihat # EL1
        u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec); 
        Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat); Z_St = t(Z_S)
        
        init = rep(0, length(Zbar)); init[length(init)] = 1
        
        nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, v_S = v_S, Z_S = Z_S, 
                              Z_St = Z_St, Zbar = Zbar, type = "EL", del = del, n = n,
                              method = "Newton", control = list(maxit = 1e5, allowSingular = T),
                              xscalm = "auto")
        if(nleqslv_res$termcd != 1){
          if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = "EL", del = del, n = n))) > 1e-5)
            w_S = NA
        }else{
          w_S = f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = "EL", del = del, n = n, returnw = T)             
        }
        
        theta_res = c(theta_res, EL0 = sum(y_S * w_S) / n) # EL1
        
        gammahat = solve(Z_St %*% (Z_S * gprime1(1 / pihat[Index_S], type = "EL", del = del) / d_S / v_S),
                         Z_St %*% (y_S * gprime1(1 / pihat[Index_S], type = "EL", del = del) / d_S / v_S))
        yhat2 = drop(cbind(1, x, u_vec) %*% gammahat)
        kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
                      t(hhat[Index_S,]) %*% ((y_S - yhat2[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        eta = yhat2 + drop(hhat %*% kappa)
        eta[Index_S] = eta[Index_S] + (y_S - yhat2[Index_S]) * w_S - drop(hhat %*% kappa)[Index_S] / pihat[Index_S]
        var_res = c(var_res, EL0 = var(eta) / n)
        
        
        d_S = rep(1, sum(delta)); #u_vec = -pihat # EL1
        v_S = vx[Index_S]; u_vec = -pihat * vx # EL2
        u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec); 
        Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat); Z_St = t(Z_S)
        
        init = rep(0, length(Zbar)); init[length(init)] = 1
        
        nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, v_S = v_S, Z_S = Z_S, 
                              Z_St = Z_St, Zbar = Zbar, type = "EL", del = del, n = n,
                              method = "Newton", control = list(maxit = 1e5, allowSingular = T),
                              xscalm = "auto")
        if(nleqslv_res$termcd != 1){
          if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = "EL", del = del, n = n))) > 1e-5)
            w_S = NA
        }else{
          w_S = f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = "EL", del = del, n = n, returnw = T)             
        }
        
        theta_res = c(theta_res, EL = sum(y_S * w_S) / n) # EL1
        
        gammahat = solve(Z_St %*% (Z_S * gprime1(1 / pihat[Index_S], type = "EL", del = del) / d_S / v_S),
                         Z_St %*% (y_S * gprime1(1 / pihat[Index_S], type = "EL", del = del) / d_S / v_S))
        yhat2 = drop(cbind(1, x, u_vec) %*% gammahat)
        kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
                      t(hhat[Index_S,]) %*% ((y_S - yhat2[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        eta = yhat2 + drop(hhat %*% kappa)
        eta[Index_S] = eta[Index_S] + (y_S - yhat2[Index_S]) * w_S - drop(hhat %*% kappa)[Index_S] / pihat[Index_S]
        var_res = c(var_res, EL = var(eta) / n)
        
        # nlmres= nlm(targetftn0, Zbar[length(Zbar)], d_S = d_S, Z_S = Z_S, Z_St = Z_St, 
        #             init = init, Zbar = Zbar, type = type, del = del)
        # # if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) print(nlmres$estimate)
        # if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) stop(nlmres$code)
        # W = nlmres$estimate
        # if(nlmres$minimum >= .Machine$double.xmax){
        #   # stop(targetftn(Zbar[length(Zbar)], d_S = d_S, Z_S = Z_S, Z_St = Z_St, init = init, Zbar = Zbar, type = type))
        #   w_S = NA
        # }else{
        #   w_S = targetftn0(W, d_S, Z_S, Z_St, init, Zbar, type, returnw = T, del = del)
        #   
        #   # nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
        #   #                       Z_St = Z_St, Zbar = c(1, colMeans(x), W), type = type, del = del,
        #   #                       method = "Newton", control = list(maxit = 1e5, allowSingular = T),
        #   #                       xscalm = "auto")
        #   # w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, returnw = T)   
        # }
        # theta_res = c(theta_res, Qin = sum(y_S * w_S) / n) # Qin
        
        # gammahat = solve(Z_St %*% (Z_S * gprime1(1 / pihat[Index_S], type, del = del) / d_S),
        #                  Z_St %*% (y_S * gprime1(1 / pihat[Index_S], type, del = del) / d_S))
        # drop(t(Zbar) %*% gammahat)
        # sum((y_S - drop(Z_S %*% gammahat)) * 1 / pihat[Index_S]) / n
        
        #### test_indices = folds[[1]];  
        type = "GE"
        targetftn = function(del, pihat, x, y, folds){
          if(is.nan(del)) return(.Machine$double.xmax)
          # print(del)
          # if(del > 0 | del < -100) return(.Machine$double.xmax)
          if(del == -1){
            type = "EL"
          }else if(del == 0){
            type = "ET"
          }else{
            type = "GE"
          }
          performance_metrics <- sapply(folds, function(test_indices) {
            train_indices <- setdiff(seq_len(n), test_indices)
            
            train_S2 = train_indices[train_indices %in% which(Index_S)]
            test_S2 = test_indices[test_indices %in% which(Index_S)]
            
            if(type == "EL"){
              u_vec0 = -pihat
            }else if(type == "ET"){
              u_vec0 = -log(pihat)
            }else{
              u_vec0 = (1 / pihat)^del / del
            }
            u_vec0 = u_vec0 * vx # EL2
            if(any(is.infinite(u_vec0))) return(.Machine$double.xmax)
            
            d_S = rep(1, length(train_S2)); # EL1
            v_S = vx[train_S2]
            
            u_vec = u_vec0[train_indices]
            u_vec_S = u_vec0[train_S2]; Uhat = mean(u_vec);
            Z_S = cbind(1, x[train_S2,], u_vec_S); Zbar = c(1, colMeans(x[train_indices,]), Uhat); Z_St = t(Z_S)
            init = rep(0, length(Zbar)); init[length(init)] = 1
            nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, v_S = v_S, 
                                  Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n,
                                  method = "Newton", control = list(maxit = 1e5, allowSingular = T),
                                  xscalm = "auto")
            # f(init, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n)
            # h(init, d_S = d_S, v_S = v_S, Z_S = Z_S, Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n)
            if(nleqslv_res$termcd != 1){
              lambda_S = nleqslv_res$x
              if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n))) > 1e-5){
                # print(lambda_S)
                return(.Machine$double.xmax)
                # stop(nleqslv_res)
              }
            }else{
              lambda_S = nleqslv_res$x
            }
            
            Amat = Z_St %*% (Z_S * gprime1(1 / pihat[train_S2], type, del = del) / d_S)
            if(rcond(Amat) < .Machine$double.eps){
              return(.Machine$double.xmax)
            }else{
              gammahat = solve(Amat,
                               Z_St %*% (y[train_S2] * gprime1(1 / pihat[train_S2], type, del = del) / d_S))
            }
            
            hmat = cbind(1, x) * pihat
            Bmat = t(hmat[train_S2,]) %*% (hmat[train_S2,] * (1 / pihat[train_S2] - 1) / pihat[train_S2])
            if(rcond(Bmat) < .Machine$double.eps){
              return(.Machine$double.xmax)
            }else{
              etahat = solve(Bmat,
                             t(hmat[train_S2,]) %*% ((
                               y[train_S2] - drop(Z_S %*% gammahat)) *
                                 (1 / pihat[train_S2] - 1) / pihat[train_S2]))
            }
            
            d_S = rep(1, length(test_S2)); # EL1
            v_S = vx[test_S2]; # EL2
            
            u_vec = u_vec0[test_indices]
            u_vec_S = u_vec0[test_S2]; Uhat = mean(u_vec);
            Z_S = cbind(1, x[test_S2,], u_vec_S); Zbar = c(1, colMeans(x[test_indices,]), Uhat); Z_St = t(Z_S)
            
            w_S = f(lambda_S, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n, returnw = T)
            # w_S[is.nan(w_S)] <- 0
            if(any(is.infinite(w_S) | is.nan(w_S))){
              # print(w_S)
              return(.Machine$double.xmax)
            }
            
            # sum((w_S * (y[test_S2] - drop(Z_S %*% gammahat)) -
            #        1 / pihat[test_S2] * (drop(hmat[test_S2,]) %*% etahat))^2)
            
            # sum((w_S * (y[test_S2] - drop(Z_S %*% gammahat)) -
            #        1 / pihat[test_S2] * (drop(hmat[test_S2,]) %*% etahat))^2 * (1 - pihat[test_S2])) +
            #   sum(w_S^2 * pihat[test_S2] * vx[test_S2])
            
            # sum(w_S^2 * (y[test_S2] - drop(Z_S %*% gammahat))^2)
            # sum(1 / pihat[test_S2] * (1 / pihat[test_S2] - 1) *
            #       (y[test_S2] - drop(Z_S%*% gammahat) - drop(hmat[test_S2,]) %*% etahat)^2)
            
            # sum(w_S * (1 / pihat[test_S2] - 1) * (y[test_S2] - drop(Z_S%*% gammahat))^2 +
            #       1 / pihat[test_S2] * (1 / pihat[test_S2] - 1) * (drop(hmat[test_S2,]) %*% etahat)^2)
            
            sum(1 / pihat[test_S2] * (1 / pihat[test_S2] - 1) * (y[test_S2] - drop(Z_S %*% gammahat) - drop(hmat[test_S2,]) %*% etahat)^2)
            # sum(w_S * (1 / pihat[test_S2] - 1) * (drop(Z_S %*% gammahat)- y[test_S2])^2)
          })
          return(ifelse(is.infinite(sum(performance_metrics)),
                        .Machine$double.xmax, sum(performance_metrics)))
        }
        
        folds <- createFolds(1:n, k = 5, list = TRUE)
        
        lower = -2; upper = 1
        xtmp = seq(from = lower, to = upper, len = (upper - lower) * 10 + 1)
        # xtmp = seq(from = lower, to = upper, len = (upper - lower) * 10)
        ytmp = sapply(xtmp, function(del)
          targetftn(del, pihat = pihat, x = x, y = y, folds = folds)
        )
        ytmp = ifelse(ytmp >= .Machine$double.xmax, Inf, ytmp)
        # plot(xtmp,log(ytmp), type = "l", main = cnt)
        
        # nlmres= nlm(targetftn, xtmp[which.min(ytmp)], pihat = pihat, x = x, y = y, folds = folds)
        # del = nlmres$estimate
        
        nlmres= nlminb(xtmp[which.min(ytmp)], targetftn,
                       pihat = pihat, x = x, y = y, folds = folds,
                       lower = lower, upper = upper)
        # nlmres$par
        (del = nlmres$par)
        
        if(del == -1){
          type = "EL"
        }else if(del == 0){
          type = "ET"
        }else{
          type = "GE"
        }
        
        if(type == "EL"){
          u_vec = -pihat
        }else if(type == "ET"){
          u_vec = -log(pihat)
        }else if(type == "GE"){
          u_vec = (1 / pihat)^del / del
        }
        u_vec = u_vec * vx # EL2
        
        # type = "GE"; del = 2 # To be removed
        d_S = rep(1, sum(delta)); # EL1
        v_S = vx[Index_S]; # EL2
        
        # type = "SL"; del = 1 ; u_vec = 1 / (pihat)# To be removed
        # type = "HD"; del = NA; u_vec = -sqrt(pihat) # To be removed
        u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec);
        Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat); Z_St = t(Z_S)
        
        init = rep(0, length(Zbar)); init[length(init)] = 1
        
        nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, v_S = v_S, 
                              Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n,
                              method = "Newton", control = list(maxit = 1e5, allowSingular = T),
                              xscalm = "auto")
        if(nleqslv_res$termcd != 1){
          if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n))) > 1e-5)
            w_S = NA
        }else{
          w_S = f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n, returnw = T)
        }
        theta_res = c(theta_res, GEC = sum(y_S * w_S) / n) # proposed
        
        gammahat = solve(Z_St %*% (Z_S * gprime1(1 / pihat[Index_S], type = "GE", del = del) / d_S / v_S),
                         Z_St %*% (y_S * gprime1(1 / pihat[Index_S], type = "GE", del = del) / d_S / v_S))
        yhat2 = drop(cbind(1, x, u_vec) %*% gammahat)
        kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
                      t(hhat[Index_S,]) %*% ((y_S - yhat2[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        eta = yhat2 + drop(hhat %*% kappa)
        eta[Index_S] = eta[Index_S] + (y_S - yhat2[Index_S]) * w_S - drop(hhat %*% kappa)[Index_S] / pihat[Index_S]
        var_res = c(var_res, GEC = var(eta) / n)
        
        CR_res = ifelse(abs(theta_res - theta) < qnorm(0.975) * sqrt(var_res), 1, 0)
        # type = "EL"; d_S = vx[Index_S]; u_vec = -pihat * vx
        # u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec);
        # Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat); Z_St = t(Z_S)
        # 
        # init = rep(0, length(Zbar)); init[length(init)] = 1
        # 
        # nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S,
        #                       Z_St = Z_St, Zbar = Zbar, type = type, del = del,
        #                       method = "Newton", control = list(maxit = 1e5, allowSingular = T),
        #                       xscalm = "auto")
        # 
        # # f(init, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del)
        # # colSums(Z_S / pihat[Index_S]) / n- Zbar
        # 
        # if(nleqslv_res$termcd != 1){
        #   if(max(abs(f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del))) > 1e-5)
        #     w_S = NA
        # }else{
        #   w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, returnw = T)
        # }
        # theta_res = c(theta_res, EL2 = sum(y_S * w_S) / n) # EL2
        
        # w = CVXR::Variable(length(d_S))
        # constraints <- list(Z_St %*% w / n == Zbar)
        # 
        # Phi_R <- Maximize(sum(log(w) * vx[Index_S]))
        # 
        # prob <- CVXR::Problem(Phi_R, constraints)
        # res <- CVXR::solve(prob, solver = "ECOS_BB")
        # 
        # w_S = drop(res$getValue(w))
        # # drop(res$getDualValue(constraints[[1]]))
        # theta_res = c(theta_res, EL2 = sum(y_S * w_S) / n) # EL2
      }
      alpha_vec = c(alpha_vec, nlmres$par)
      theta_mat = cbind(theta_mat, theta_res - theta)
      var_mat = cbind(var_mat, var_res)
      CR_mat = cbind(CR_mat, CR_res)
    }
    list(theta_mat, alpha = alpha_vec, var_mat, CR_mat)
  }
# final_res
if(!interactive()) save.image(paste(timenow0, ".RData", sep = ""))
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

res1 <- matrix(rowMeans(do.call(cbind, lapply(final_res1, c)), na.rm = TRUE), 
               ncol = ncol(final_res1[[1]]))
colnames(res1) = rnames; row.names(res1) = row.names(final_res1[[1]])

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

res1
res2

xtable::xtable(res1, digits = 3)
xtable::xtable(res2, digits = 3)

Ncol = ncol(final_res2[[1]])
Nrow = nrow(final_res2[[1]])

for(ncols in 1:Ncol){
  boxres = t(do.call(cbind, lapply(final_res1, c)))[, ((ncols - 1) * Nrow + 1) : (ncols * Nrow)] + theta
  colnames(boxres) <- row.names(final_res1[[1]])
  if(!interactive()) png(paste("boxplot", ncols, ".png", sep = ""))
  boxplot(boxres, main = paste(paste(c("RP", "OR"), ifelse(modelcases[ncols,], "C", "M"), sep = ":"), collapse = ", "))
  abline(h = theta, col = "red")
  if(!interactive()) dev.off()
}


if(!interactive()) png("boxplots.png", width = 960, height = 960)
op <- par(); par(mfrow = c(2,2), mar = rep(3, 4))
for(ncols in 1:Ncol){
  boxres = t(do.call(cbind, lapply(final_res1, c)))[, ((ncols - 1) * Nrow + 1) : (ncols * Nrow)] + theta
  colnames(boxres) <- row.names(final_res1[[1]])
  boxplot(boxres, main = paste(paste(c("RP", "OR"), ifelse(modelcases[ncols,], "C", "M"), sep = ":"), collapse = ", "), 
          ylim = c(theta - 10, theta + 10), cex.axis=1.2, cex.main = 1.2)
  # boxplot(boxres, main = rnames[ncols], ylim = c(theta - 10, theta + 10), cex.axis=0.8)
  abline(h = theta, col = "red")
}
par(op)
if(!interactive()) dev.off()

final_res_alpha = lapply(final_res, function(x) x[[2]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
alpha_df <- as.data.frame(do.call("rbind", final_res_alpha))
colnames(alpha_df) <- rnames
alpha_df_long <- pivot_longer(alpha_df, cols = everything(), names_to = "Column", values_to = "Value")
alpha_df_long$Column <- factor(alpha_df_long$Column, levels = rnames)

# Plot using ggplot2
pGG <- ggplot(alpha_df_long, aes(x = Value, fill = Column)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "Pastel1") + # Optional: use a color palette
  theme_minimal() + # Optional: use a minimal theme for the plot
  labs(title = "Density plot of alpha", x = "alpha", y = "Density") +
  theme(legend.title = element_blank()) # Optional: remove legend title
pGG

if(!interactive()) ggsave("hist_alpha.png", plot = pGG)


gsub("\\\\addlinespace", "", kable(res1, "latex", booktabs = TRUE, digits = 3) %>%
       kable_styling())

gsub("\\\\addlinespace", "", kable(res2, "latex", booktabs = TRUE, digits = 3) %>%
       kable_styling())

final_res_var = lapply(final_res, function(x) x[[3]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
res_var <- matrix(rowMeans(do.call(cbind, lapply(final_res_var, c)), na.rm = TRUE), 
                  ncol = ncol(final_res_var[[3]]))
colnames(res_var) = rnames; row.names(res_var) = row.names(final_res_var[[3]])

(res_var - res2^2) / res2^2

final_res_CR = lapply(final_res, function(x) x[[4]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
res_CR <- matrix(rowMeans(do.call(cbind, lapply(final_res_CR, c)), na.rm = TRUE), 
                 ncol = ncol(final_res_CR[[4]]))
colnames(res_CR) = rnames; row.names(res_CR) = row.names(final_res_CR[[4]])
# res_CR

gsub("\\\\addlinespace", "", kable((res_var - res2^2) / res2^2, "latex", booktabs = TRUE, digits = 3) %>%
       kable_styling())

gsub("\\\\addlinespace", "", kable(res_CR, "latex", booktabs = TRUE, digits = 3) %>%
       kable_styling())
