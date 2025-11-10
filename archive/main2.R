# High-dimensional Simulation setup

if(!interactive()){
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
}else{
  args <- c(300)
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# install.packages( "MatrixModels", type="win.binary" )  

library(nleqslv)
library(CVXR)
library(RCAL)
library(CBPS)
library(statmod)
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
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
registerDoRNG(seed = 11)

# modelcases = expand.grid(c(T,F),c(T,F),c(T,F))
modelcases = expand.grid(c(T,F),c(T,F))
rnames <- apply(apply(modelcases, 2, ifelse, "C", "M"), 1, paste, collapse = "")

# theta = 210 # Kang & Schaffer

final_res <- foreach(
  simnum = 1:SIMNUM, 
  .packages = c("nleqslv", "CVXR", "caret", "RCAL", "CBPS", "glmnet"),
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
      # Wang & Kim 2024, modified
      n = 500 # n = 200
      pcol = 100; # pcol = 4 is True
      # rho <- 0.3 # AR(1) correlation parameter
      # L <- chol(outer(1:pcol, 1:pcol, function(j, k) rho^abs(j - k)))
      # x <- matrix(rnorm(n * pcol, 2, 1), nrow = n, ncol = pcol) %*% L
      x = matrix(rnorm(n * pcol, 2, 1), nc= pcol)
      # x = ifelse(x < 0, 0, x)
      # x = ifelse(x > 4, 4, x)
      z = x
      if(RM  == T){
        pi = 1 / (1 + exp(-(-.5 - x[,3] + 0.5 * x[,4] + 0.5 * x[,5] - 0.25 * x[,6])))
      }else if(RM == F){
        #   pi = 1 / (1 + exp(-(- 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.75 * (x[,2] - 2)^2 -
        #                         0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
        # pi = ifelse(pi < 0.05, 0.05, pi)
        # pi = ifelse(pi > 0.7, 0.7, pi)
          # pi = 1 / (1 + exp(-(1 - 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
          #                       0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
          # pi = pt(- 0.1 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
          #           0.2 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2, 3)
        pi = 1 / (1 + exp(-(-.5 - 0.5 * (x[,3] - 2) * (x[,4] - 5) - 0.5 * (x[,4] - 2.5)^2 -
                              0.5 * (x[,5] - 1.5)^2 * (x[,6] - 2)^2 )))
          pi = ifelse(pi < 0.05, 0.05, pi)
          # pi = ifelse(pi > 0.7, 0.7, pi)
          summary(pi)
          mean(pi)
      }
      vx = rep(1, n) # original
      e = rnorm(n, 0, sqrt(vx))
      
      if(OM == T){
        y = 1 + x[,1] - x[,2] + x[,3] - x[,4] + e
        theta = 1
      }else if(OM == F){
        # y = 1 - 0.5 * (x[,1] - 2)^2 * (x[,2]-5)^2 + 0.2 * (x[,2] - 1)^2 * (x[,3] - 1.5)^3 -
        #   0.5 * cos(x[,3]^2) * sin(x[,4]^2) + e
        # theta = -3.351918
          # y = 1 + 0.5 * x[,1]^2 * x[,2]^2 + 0.2 * x[,2]^2 * x[,3]^3 +
          #   0.5 * cos(x[,3]^2) + 0.5 * sin(x[,4]^2) + e
          # theta = 27.58047  
          # y = 1 - 0.5 * (x[,1] - 3)^2 * (x[,2]-5)^2 + 0.2 * (x[,2] - 1)^2 * (x[,3] - 1.5)^3 -
          #   0.5 * cos(x[,3]^2) - 0.5 * sin(x[,4]^2) + e
          # theta = -8.430471
        
        # gh <- gauss.quad(100, "hermite")
        # nodes <- sqrt(2) * gh$nodes + 2
        # weights <- gh$weights / sqrt(base::pi)
        # 
        # grid <- setNames(expand.grid(nodes, nodes), c("x1", "x2"))
        # w_grid <- expand.grid(weights, weights)
        # f1 <- function(x1, x2) 1 - (x1- 2) * (x2 - 1) + 0.75 * (x1 - 1.5)^2 -sin(x2^2)
        # sum(do.call(f1, grid) * Reduce(`*`, w_grid))
        
        y = 1 - (x[,1] - 2) * (x[,2] - 1) + 0.75 * (x[,3] - 1.5)^2 - sin(x[,4]^2)  + e
        theta = 1.805657  # OM1
      }
      
      theta_res = NULL
      var_res = NULL
      CR_res = NULL
      
      delta = rbinom(n, 1, pi)
      Index_S = (delta == 1)
      
      y_S = y[Index_S]
      x_S = x[Index_S,]
      z_S = z[Index_S,]
      # x1 = cbind(1, x)
      
      del = NA
      type = "EL"
      
      data = data.frame(y, delta, vx = vx, x = x, z = z)
      data_S = data[Index_S,]
      
      col_RM = 1:pcol    
      xR_S = x_S[, col_RM]; xR = x[,col_RM]
      
      # library(RCAL)
      # mn.cv.RCAL <- mn.regu.cv(fold=5*c(1,1), nrho=(1+10)*c(1,1), rho.seq=NULL, y, delta, x, 
      #            ploss="cal", yloss="gaus")
      # pihat.RCAL <- mn.cv.RCAL$ps$sel.fit[,1] # pihat using RCAL
      # mn.cv.RCAL$ps$sel.fit[,2] # yhat using RCAL
      # names(mn.cv.RCAL$ps)
      # mn.cv.RCAL$ps$sel.rho
      # mn.cv.RCAL$or$sel.rho
      # mn.cv.RCAL$est$ipw # IPW est
      # mn.cv.RCAL$est$est # AIPW est

      ps.cv.RCAL <- glm.regu.cv(fold=3, nrho=1+3, y=delta, x=xR, loss="cal")
      pihat.RCAL <- ps.cv.RCAL$sel.fit[,1]
      
      # ps.cv.RCAL <- glm.regu.cv(fold=5, nrho=1+10, y=delta, x=xR, loss="cal")
      # ps.cv.RCAL$sel.fit[,1] # pihat.cv.RCAL
      # ps.RCAL = glm.regu(delta, x1, rhos = rep(ps.cv.RCAL$sel.rho[1], ncol(x1)), loss = "cal") # Non default pihat.RCAL
      # ps.RCAL = glm.regu(delta, xR, rhos = rep(ps.cv.RCAL$sel.rho[2], ncol(xR)), loss = "cal") # ps.cv.RCAL$sel.rho[2] is default
      # ps.RCAL = glm.regu(delta, xR, rhos = rep(0.1, ncol(xR)), loss = "cal")
            
      # phihat.RCAL = c(ps.RCAL$inter, ps.RCAL$bet)
      
      # pihat.RCAL <- ps.RCAL$fit # pihat.cral
      # pihat.RCAL <- 1 / findphi0(phihat.RCAL, xR = xR, Z = xR, delta = delta, returnw = T)

      # sum(y_S / pihat.RCAL[Index_S]) / sum(1 / pihat.RCAL[Index_S]) #IPW.RCAL
      
      ps.cv.MLE <- glm.regu.cv(fold=3, nrho=1+3, y=delta, x=xR, loss="ml")
      pihat.MLE <- ps.cv.MLE$sel.fit[,1] # pihat using MLE
      # ps.MLE = glm.regu(delta, xR, rhos = rep(ps.cv.MLE$sel.rho[2], ncol(xR)), loss = "ml")
      
      # ps.MLE = glm.regu(delta, xR, rhos = rep(0.02, ncol(xR)), loss = "ml")
      # phihat.MLE = c(ps.MLE$inter, ps.MLE$bet)
      
      # pihat.MLE <- ps.MLE$fit # pihat.mle

      # library(CBPS)
      # ps.CBPS <- hdCBPS(reformulate(paste0("x.", col_RM), response = "delta"), data = data,
      #        y = y, ATT = 0)
      # pihat.CBPS <- ps.CBPS$fitted.values # pihat.CBPS
      # 
      # sum(y_S / pihat.CBPS[Index_S]) / sum(1 / pihat.CBPS[Index_S]) #IPW.cbps

      # library(glmnet)

      col_OM = 1:pcol
      xO_S = x_S[, col_OM]; xO = x[, col_OM]
      # Omodel.cv = cv.glmnet(xO_S, y_S, weights = 1 / vx[Index_S])
      # betahat.glmnet = as.vector(coef(Omodel.cv, s = "lambda.min"))
      # Omodel.cv$glmnet.fit
      # yhat = drop(predict(Omodel.cv, newx = xO, s = "lambda.min"))
      
      # plot(drop(predict(Omodel.cv, newx = xO, s = "lambda.min")), 1 + x[,1] - x[,2] + x[,3] - x[,4])

      or.cv <- glm.regu.cv(fold=3, nrho=1+3, y = y_S, x = xO_S, loss="gaus")
      yhat <- c( cbind(1,x)%*%or.cv$sel.bet[,1] )
      # if(OM){
      #   sel.rho = 0.6
      # }else if(OM){
      #   sel.rho = 0.2
      # } 
      # # or <- glm.regu(y = y_S, x = xO_S, rhos = rep(or.cv$sel.rho[1], ncol(xR)), loss="gaus")
      # or <- glm.regu(y = y_S, x = cbind(1, xO_S), rhos = rep(0.2, ncol(xR)+1), loss="gaus")
      # yhat <- c( cbind(1,xO) %*% or$bet )
      
      
      
      
      for(pimethod in c(1,3)){
        if(pimethod == 1){
          pihat = pihat.MLE
          hhat = pihat * xR # MLE
        }else if(pimethod == 2){
          pihat = pihat.CBPS
          hhat = 1 / (1 - pihat) * xR
        }else if(pimethod == 3){
          pihat = pihat.RCAL
          # hhat = pihat / (1 - pihat) * xR # Tan
          hhat = xR # Tan
        } 
        mn.aipw.res = mn.aipw(y, delta, fp=pihat, fo=yhat)
        
        # y_IPW = sum(y_S / pihat[Index_S]) / sum(1 / pihat[Index_S])
        y_IPW = mn.aipw.res$ipw
        theta_res = c(theta_res, IPW = y_IPW)  # IPW
        # if (rcond(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S])) < 1e-12) {
        #   kappa = rep(NA, ncol(hhat))
        # } else {
        #   kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
        #                 t(hhat[Index_S,]) %*% ((y_S - y_IPW) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        # }
        # 
        # eta = drop(hhat %*% kappa)
        # eta[Index_S] = eta[Index_S] + (y_S - y_IPW - eta[Index_S]) / pihat[Index_S]
        # var_res = c(var_res, IPW = var(eta) / n)
        var_res = c(var_res, IPW = NA)  

        y_AIPW = mn.aipw.res$est
        # y_AIPW = (sum(yhat) / n + sum((y_S - yhat[Index_S]) / pihat[Index_S]) / n)
        # y_AIPW = (sum(yhat) / n + sum((y_S - yhat[Index_S]) / pihat[Index_S]) / sum(1 / pihat[Index_S]))
        
        theta_res = c(theta_res, AIPW = y_AIPW) # AIPW
        # if (rcond(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S])) < 1e-12) {
        #   gammahat = rep(NA, ncol(hhat[Index_S,]))
        # } else {
        #   kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
        #                 t(hhat[Index_S,]) %*% ((y_S - yhat[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        # }
        # 
        # eta = yhat + drop(hhat %*% kappa)
        # eta[Index_S] = eta[Index_S] + (y_S - eta[Index_S]) / pihat[Index_S]
        var_res = c(var_res, AIPW = mn.aipw.res$var)
        
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
        
        for(type in c("EL", "ET", "GE")){
          print(type)
          if(type == "GE") del = -0.5
          d_S = rep(1, sum(delta)); 
          v_S = rep(1, sum(delta)); 
          if(type == "EL"){
            u_vec_tmp = -pihat
          }else if(type == "ET"){
            u_vec_tmp = -log(pihat)
          }else if(type == "GE"){
            u_vec_tmp = (1 / pihat)^del / del
          }
          u_vec = matrix(u_vec_tmp, ncol = 1) # EL1
          u_vec_S = u_vec[Index_S,]; Uhat = colMeans(u_vec); 
          Z_S = cbind(1, u_vec_S, yhat[Index_S])
          Zbar = c(1, Uhat, mean(yhat)); Z_St = t(Z_S)
          
          # w = CVXR::Variable(length(d_S))
          # constraints <- list(Z_St %*% w / n == Zbar, 
          #                     abs(t(xO_S) %*% w / n - colMeans(xO)) <= 1)
          # 
          # # tail(sort(abs(t(xO_S) %*% (1 / pihat[Index_S]) / n - colMeans(xO))), 300)
          # 
          # if(type == "EL"){
          #   Phi_R <- Maximize(sum(log(w)))
          # }else if(type == "ET"){
          #   Phi_R <- Maximize(sum(entr(w) - w))
          # }else if(type == "GE"){
          #   Phi_R <- Maximize(sum(-power(w, del + 1) * del^(-1) * (del + 1)^(-1)))          
          # }
          # 
          # prob <- CVXR::Problem(Phi_R, constraints)
          # (res <- CVXR::solve(prob, solver = "ECOS_BB"))
          # if(res$status == "solver_error"){
          #   w_S = rep(NA, length(d_S))
          # }else{
          #   w_S = drop(res$getValue(w))
          # }
          # #  max(abs(Z_St %*% w_S / n - Zbar)); max(abs(t(xO_S) %*% w_S / n - colMeans(xO)))
          # 
          # theta_res = c(theta_res, setNames(sum(y_S * w_S) / n, paste(type, 0, sep = ""))) # EL1
          # var_res = c(var_res, setNames(NA, paste(type, 0, sep = "")))
          

          
          
          d_S = rep(1, sum(delta)); #u_vec = -pihat # EL1
          vx = rep(1, length(delta)); v_S = vx[Index_S];
          u_vec = matrix(u_vec_tmp * vx, ncol = 1) # EL1
          u_vec_S = u_vec[Index_S,]; Uhat = colMeans(u_vec); 
          Z_S = cbind(1, u_vec_S, yhat[Index_S]); Zbar = c(1, Uhat, mean(yhat)); Z_St = t(Z_S)          
          Z = cbind(1, u_vec, yhat)
          orth_vec = -xR * ((1 - pihat) / pihat / gprime1(1 / pihat, type = type, del = del) * vx); orth_vec_S = orth_vec[Index_S,]
          
          w = CVXR::Variable(length(d_S))
          constraints <- list(Z_St %*% w / n == Zbar)
          
          if(type == "EL"){
            Phi_R <- Maximize(sum(log(w)))
          }else if(type == "ET"){
            Phi_R <- Maximize(sum(entr(w) - w))
          }else if(type == "GE"){
            Phi_R <- Maximize(sum(-power(w, del + 1) * del^(-1) * (del + 1)^(-1)))          
          }
          
          prob <- CVXR::Problem(Phi_R, constraints)
          res <- CVXR::solve(prob, solver = "ECOS_BB")
          if(res$status == "solver_error"){
            w_S = rep(NA, length(d_S))
          }else{
            w_S = drop(res$getValue(w))
          }
          
          theta_res = c(theta_res, setNames(sum(y_S * w_S) / n, paste(type, 0, sep = ""))) # EL1
          var_res = c(var_res, setNames(NA, paste(type, 0, sep = "")))
          # var_res = c(var_res, setNames(sum((c(y_S * w_S + yhat[Index_S] * (1 - w_S), 
          #                                      yhat[!Index_S]) - sum(y_S * w_S) / n)^2) / n / (n-1), paste(type, 0, sep = "")))          
          
          # max(abs(Z_St %*% w_S / n - Zbar));
          tau1 = max(abs(t(xR_S) %*% w_S / n - colMeans(xR))); 
          tau2 = max(abs(t(orth_vec_S) %*% w_S / n - colMeans(orth_vec))); 
          
          constraints <- list(Z_St %*% w / n == Zbar, 
                              abs(t(xR_S) %*% w / n - colMeans(xR)) <= 0.7 * tau1,
                              abs(t(orth_vec_S) %*% w / n - colMeans(orth_vec)) <= 0.7 * tau2)
          
          if(type == "EL"){
            Phi_R <- Maximize(sum(log(w)))
          }else if(type == "ET"){
            Phi_R <- Maximize(sum(entr(w) - w))
          }else if(type == "GE"){
            Phi_R <- Maximize(sum(-power(w, del + 1) * del^(-1) * (del + 1)^(-1)))          
          }
          
          prob <- CVXR::Problem(Phi_R, constraints)
          res <- CVXR::solve(prob, solver = "ECOS_BB")
          
          if(res$status == "solver_error"){
            w_S = rep(NA, length(d_S))
          }else{
            w_S = drop(res$getValue(w))
          }
          
          # max(abs(t(orth_vec_S) / pihat[Index_S] / n - colMeans(orth_vec)))
          
          # max(abs(t(xR_S) / pihat[Index_S] / n - colMeans(xR)))
          
           # max(abs(Z_St %*% w_S / n - Zbar));max(abs(t(xR_S) %*% w_S / n - colMeans(xR))); max(abs(t(orth_vec_S) %*% w_S / n - colMeans(orth_vec))); 
          theta_res = c(theta_res, setNames(sum(y_S * w_S) / n, type)) # EL1
          
          var_res = c(var_res, setNames(sum((c(y_S * w_S + yhat[Index_S] * (1 - w_S), 
                                               yhat[!Index_S]) - sum(y_S * w_S) / n)^2) / n / (n-1), type))
        }
        

        
        # kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
        #               t(hhat[Index_S,]) %*% ((y_S - yhat2[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        # eta = yhat2 + drop(hhat %*% kappa)
        # eta[Index_S] = eta[Index_S] + (y_S - yhat2[Index_S]) * w_S - drop(hhat %*% kappa)[Index_S] / pihat[Index_S]
        # var_res = c(var_res, EL = var(eta) / n)
        
        # #### test_indices = folds[[1]];  
        # type = "GE"
        # targetftn = function(del, pihat, x, y, folds){
        #   if(is.nan(del)) return(.Machine$double.xmax)
        #   # print(del)
        #   # if(del > 0 | del < -100) return(.Machine$double.xmax)
        #   if(del == -1){
        #     type = "EL"
        #   }else if(del == 0){
        #     type = "ET"
        #   }else{
        #     type = "GE"
        #   }
        #   performance_metrics <- sapply(folds, function(test_indices) {
        #     train_indices <- setdiff(seq_len(n), test_indices)
        #     
        #     train_S2 = train_indices[train_indices %in% which(Index_S)]
        #     test_S2 = test_indices[test_indices %in% which(Index_S)]
        #     
        #     if(type == "EL"){
        #       u_vec_tmp = -pihat
        #     }else if(type == "ET"){
        #       u_vec_tmp = -log(pihat)
        #     }else{
        #       u_vec_tmp = (1 / pihat)^del / del
        #     }
        #     u_vec0 = matrix(u_vec_tmp * vx, ncol = 1) # EL1
        #     u_vec0 = cbind(-xO * ((1 - pihat) / pihat / gprime1(1 / pihat, type = type, del = del) * vx), u_vec0)
        #     if(any(is.infinite(u_vec0))) return(.Machine$double.xmax)
        #     
        #     d_S = rep(1, length(train_S2)); # EL1
        #     v_S = vx[train_S2]
        #     
        #     u_vec = u_vec0[train_indices,]
        #     u_vec_S = u_vec0[train_S2,]; Uhat = colMeans(u_vec);
        #     Z_S = cbind(1, xO[train_S2,], u_vec_S); Zbar = c(1, colMeans(xO[train_indices,]), Uhat); Z_St = t(Z_S)
        #     init = rep(0, length(Zbar)); init[length(init)] = 1
        #     nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, v_S = v_S, 
        #                           Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n,
        #                           method = "Newton", control = list(maxit = 1e5, allowSingular = T),
        #                           xscalm = "auto")
        #     # f(init, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n)
        #     # h(init, d_S = d_S, v_S = v_S, Z_S = Z_S, Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n)
        #     if(nleqslv_res$termcd != 1){
        #       lambda_S = nleqslv_res$x
        #       if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n))) > 1e-5){
        #         # print(lambda_S)
        #         return(.Machine$double.xmax)
        #         # stop(nleqslv_res)
        #       }
        #     }else{
        #       lambda_S = nleqslv_res$x
        #     }
        #     
        #     Amat = Z_St %*% (Z_S * gprime1(1 / pihat[train_S2], type, del = del) / d_S)
        #     if(rcond(Amat) < .Machine$double.eps){
        #       return(.Machine$double.xmax)
        #     }else{
        #       gammahat = solve(Amat,
        #                        Z_St %*% (y[train_S2] * gprime1(1 / pihat[train_S2], type, del = del) / d_S))
        #     }
        #     
        #     hmat = cbind(1, xO) * pihat
        #     Bmat = t(hmat[train_S2,]) %*% (hmat[train_S2,] * (1 / pihat[train_S2] - 1) / pihat[train_S2])
        #     if(rcond(Bmat) < .Machine$double.eps){
        #       return(.Machine$double.xmax)
        #     }else{
        #       etahat = solve(Bmat,
        #                      t(hmat[train_S2,]) %*% ((
        #                        y[train_S2] - drop(Z_S %*% gammahat)) *
        #                          (1 / pihat[train_S2] - 1) / pihat[train_S2]))
        #     }
        #     
        #     d_S = rep(1, length(test_S2)); # EL1
        #     v_S = vx[test_S2]; # EL2
        #     
        #     u_vec = u_vec0[test_indices,]
        #     u_vec_S = u_vec0[test_S2,]; Uhat = colMeans(u_vec);
        #     Z_S = cbind(1, xO[test_S2,], u_vec_S); Zbar = c(1, colMeans(xO[test_indices,]), Uhat); Z_St = t(Z_S)
        #     
        #     w_S = f(lambda_S, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n, returnw = T)
        #     # w_S[is.nan(w_S)] <- 0
        #     if(any(is.infinite(w_S) | is.nan(w_S))){
        #       # print(w_S)
        #       return(.Machine$double.xmax)
        #     }
        #     
        #     # sum((w_S * (y[test_S2] - drop(Z_S %*% gammahat)) -
        #     #        1 / pihat[test_S2] * (drop(hmat[test_S2,]) %*% etahat))^2)
        #     
        #     # sum((w_S * (y[test_S2] - drop(Z_S %*% gammahat)) -
        #     #        1 / pihat[test_S2] * (drop(hmat[test_S2,]) %*% etahat))^2 * (1 - pihat[test_S2])) +
        #     #   sum(w_S^2 * pihat[test_S2] * vx[test_S2])
        #     
        #     # sum(w_S^2 * (y[test_S2] - drop(Z_S %*% gammahat))^2)
        #     # sum(1 / pihat[test_S2] * (1 / pihat[test_S2] - 1) *
        #     #       (y[test_S2] - drop(Z_S%*% gammahat) - drop(hmat[test_S2,]) %*% etahat)^2)
        #     
        #     # sum(w_S * (1 / pihat[test_S2] - 1) * (y[test_S2] - drop(Z_S%*% gammahat))^2 +
        #     #       1 / pihat[test_S2] * (1 / pihat[test_S2] - 1) * (drop(hmat[test_S2,]) %*% etahat)^2)
        #     
        #     sum(1 / pihat[test_S2] * (1 / pihat[test_S2] - 1) * (y[test_S2] - drop(Z_S %*% gammahat) - drop(hmat[test_S2,]) %*% etahat)^2)
        #     # sum(w_S * (1 / pihat[test_S2] - 1) * (drop(Z_S %*% gammahat)- y[test_S2])^2)
        #   })
        #   return(ifelse(is.infinite(sum(performance_metrics)),
        #                 .Machine$double.xmax, sum(performance_metrics)))
        # }
        # 
        # folds <- createFolds(1:n, k = 5, list = TRUE)
        # 
        # lower = -2; upper = 1
        # xtmp = seq(from = lower, to = upper, len = (upper - lower) * 10 + 1)
        # # xtmp = seq(from = lower, to = upper, len = (upper - lower) * 10)
        # ytmp = sapply(xtmp, function(del)
        #   targetftn(del, pihat = pihat, x = xO, y = y, folds = folds)
        # )
        # ytmp = ifelse(ytmp >= .Machine$double.xmax, Inf, ytmp)
        # # plot(xtmp,log(ytmp), type = "l", main = cnt)
        # 
        # # nlmres= nlm(targetftn, xtmp[which.min(ytmp)], pihat = pihat, x = xO, y = y, folds = folds)
        # # del = nlmres$estimate
        # 
        # nlmres= nlminb(xtmp[which.min(ytmp)], targetftn,
        #                pihat = pihat, x = xO, y = y, folds = folds,
        #                lower = lower, upper = upper)
        # # nlmres$par
        # (del = nlmres$par)
        # 
        # if(del == -1){
        #   type = "EL"
        # }else if(del == 0){
        #   type = "ET"
        # }else{
        #   type = "GE"
        # }
        # 
        # if(type == "EL"){
        #   u_vec_tmp = -pihat
        # }else if(type == "ET"){
        #   u_vec_tmp = -log(pihat)
        # }else if(type == "GE"){
        #   u_vec_tmp = (1 / pihat)^del / del
        # }
        # u_vec = matrix(u_vec_tmp * vx, ncol = 1) # EL1
        # u_vec = cbind(-xR * ((1 - pihat) / pihat / gprime1(1 / pihat, type = type, del = del) * vx), u_vec)
        # 
        # # type = "GE"; del = 2 # To be removed
        # d_S = rep(1, sum(delta)); # EL1
        # v_S = vx[Index_S]; # EL2
        # 
        # # type = "SL"; del = 1 ; u_vec = 1 / (pihat)# To be removed
        # # type = "HD"; del = NA; u_vec = -sqrt(pihat) # To be removed
        # 
        # u_vec_S = u_vec[Index_S,]; Uhat = colMeans(u_vec);
        # Z_S = cbind(1, xO_S, u_vec_S); Zbar = c(1, colMeans(xO), Uhat); Z_St = t(Z_S)
        # 
        # init = rep(0, length(Zbar)); init[length(init)] = 1
        # 
        # nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, v_S = v_S, 
        #                       Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n,
        #                       method = "Newton", control = list(maxit = 1e5, allowSingular = T),
        #                       xscalm = "auto")
        # if(nleqslv_res$termcd != 1){
        #   if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n))) > 1e-5)
        #     w_S = NA
        # }else{
        #   w_S = f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n, returnw = T)
        # }
        # theta_res = c(theta_res, GEC = sum(y_S * w_S) / n) # proposed
        # 
        # gammahat = solve(Z_St %*% (Z_S * gprime1(1 / pihat[Index_S], type = "GE", del = del) / d_S / v_S),
        #                  Z_St %*% (y_S * gprime1(1 / pihat[Index_S], type = "GE", del = del) / d_S / v_S))
        # yhat2 = drop(cbind(1, xO, u_vec) %*% gammahat)
        # 
        # var_res = c(var_res, GEC = sum((c(y_S * w_S + drop(Z_S %*% gammahat) * (1 - w_S), 
        #                                   yhat2[!Index_S]) - sum(y_S * w_S) / n)^2) / n / (n-1))
        # # kappa = solve(t(hhat[Index_S,]) %*% (hhat[Index_S,] * (1 / pihat[Index_S] - 1) / pihat[Index_S]),
        # #               t(hhat[Index_S,]) %*% ((y_S - yhat2[Index_S]) * (1 / pihat[Index_S] - 1) / pihat[Index_S]))
        # # eta = yhat2 + drop(hhat %*% kappa)
        # # eta[Index_S] = eta[Index_S] + (y_S - yhat2[Index_S]) * w_S - drop(hhat %*% kappa)[Index_S] / pihat[Index_S]
        # # var_res = c(var_res, GEC = var(eta) / n)
        
        CR_res = ifelse(abs(theta_res - theta) < qnorm(0.975) * sqrt(var_res), 1, 0)
        # type = "EL"; d_S = vx[Index_S]; u_vec = -pihat * vx
        # u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec);
        # Z_S = cbind(1, xO_S, u_vec_S); Zbar = c(1, colMeans(xO), Uhat); Z_St = t(Z_S)
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
      nlmres = NULL
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

# xtable::xtable(res1, digits = 3)
# xtable::xtable(res2, digits = 3)

Ncol = ncol(final_res2[[1]])
Nrow = nrow(final_res2[[1]])

# for(ncols in 1:Ncol){
#   boxres = t(do.call(cbind, lapply(final_res1, c)))[, ((ncols - 1) * Nrow + 1) : (ncols * Nrow)] + theta
#   colnames(boxres) <- row.names(final_res1[[1]])
#   if(!interactive()) png(paste("boxplot", ncols, ".png", sep = ""))
#   boxplot(boxres, main = paste(paste(c("RP", "OR"), ifelse(modelcases[ncols,], "C", "M"), sep = ":"), collapse = ", "))
#   abline(h = theta, col = "red")
#   if(!interactive()) dev.off()
# }


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

# final_res_alpha = lapply(final_res, function(x) x[[2]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
# alpha_df <- as.data.frame(do.call("rbind", final_res_alpha))
# colnames(alpha_df) <- rnames
# alpha_df_long <- pivot_longer(alpha_df, cols = everything(), names_to = "Column", values_to = "Value")
# alpha_df_long$Column <- factor(alpha_df_long$Column, levels = rnames)
# 
# # Plot using ggplot2
# pGG <- ggplot(alpha_df_long, aes(x = Value, fill = Column)) +
#   geom_density(alpha = 0.5) +
#   scale_fill_brewer(palette = "Pastel1") + # Optional: use a color palette
#   theme_minimal() + # Optional: use a minimal theme for the plot
#   labs(title = "Density plot of alpha", x = "alpha", y = "Density") +
#   theme(legend.title = element_blank()) # Optional: remove legend title
# pGG
# 
# if(!interactive()) ggsave("hist_alpha.png", plot = pGG)


gsub("\\\\addlinespace", "", kable(cbind(res1, res2), "latex", booktabs = TRUE, digits = 3) %>%
  kable_styling())

final_res_var = lapply(final_res, function(x) x[[3]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
res_var <- matrix(rowMeans(do.call(cbind, lapply(final_res_var, c)), na.rm = TRUE), 
               ncol = ncol(final_res_var[[3]]))
colnames(res_var) = rnames; row.names(res_var) = row.names(final_res_var[[3]])

na_res = matrix(colSums(apply(do.call(cbind, lapply(final_res_var, c)), 1, is.na)),
                ncol = ncol(final_res_var[[1]]))
colnames(na_res) = rnames; row.names(na_res) = row.names(final_res1[[1]])
na_res

(res_var - res3^2) / res3^2

final_res_CR = lapply(final_res, function(x) x[[4]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
res_CR <- matrix(rowMeans(do.call(cbind, lapply(final_res_CR, c)), na.rm = TRUE), 
                  ncol = ncol(final_res_CR[[4]]))
colnames(res_CR) = rnames; row.names(res_CR) = row.names(final_res_CR[[4]])
# res_CR

gsub("\\\\addlinespace", "", kable(cbind((res_var - res3^2) / res3^2, res_CR), "latex", booktabs = TRUE, digits = 3) %>%
       kable_styling())

