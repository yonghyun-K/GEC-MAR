# High-dimensional Simulation setup

if(!interactive()){
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
}else{
  args <- c(100)
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

load("main2_nhanes2.Rdata")
# install.packages( "MatrixModels", type="win.binary" )  

library(nleqslv)
library(CVXR)
library(RCAL)
library(CBPS)
library(glmnet)
library(statmod)
library(dplyr)
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
library(caret)
library(tidyverse)
library(kableExtra)
library(Rmosek)

set.seed(11)
SIMNUM = args[1]

# n = nrow(nhanes0_comp)
# Index_tmp = sample(1:n, 500)
# nhanes0_comp = nhanes0_comp[Index_tmp,]
# pi = pi[Index_tmp]; delta = delta[Index_tmp]

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
# modelcases = expand.grid(c(T,F),c(T,F))
modelcases = expand.grid(c(T,F))
rnames <- apply(apply(modelcases, 2, ifelse, "C", "M"), 1, paste, collapse = "")
rnames = T
# theta = 210 # Kang & Schaffer

delta <- as.numeric(is.na(nhanes0$BPXSY1))
nhanes0_comp <- (nhanes0_comp[,c(7:32)] %>% dplyr::select(!c(RIDRETH3, DMDEDUC2, DMDHHSIZ, INDHHIN2, ALQ121)))
# nhanes0_comp <- (nhanes0_comp %>% select(!c(URXCRS, ALQ151, PEASCCT1)))
# Xtmp <- model.matrix(~ ., nhanes0_comp %>% select(!c(BPXSY1)))[, -1]
Xtmp <- model.matrix(~ ., nhanes0_comp %>% dplyr::select(!c(BPXSY1, BPXSY2, BPXSY3, BPXDI1, BPXDI2, BPXDI3)))[, -1]
cvfit <- cv.glmnet(Xtmp, delta, alpha=1, family="binomial", nfolds=10)
b <- drop(coef(cvfit, s=.03))[-1]; 
v <- names(b[b!=0]); length(v)
model1 <- glm(delta~., data=data.frame(delta, Xtmp[,v,drop=FALSE]), family=binomial)
pi <- predict(model1, type="response")

Z
# rho_vec = seq(from = 0.5, to = 0.1, by = -0.4)
# rho_vec = c(0.75, 0.5)
# rho_vec = 0.5
theta = mean(nhanes0_comp$BPXSY1)

final_res <- foreach(
  simnum = 1:SIMNUM, 
  .packages = c("nleqslv", "CVXR", "caret", "RCAL", "CBPS", "glmnet", "dplyr", "Rmosek"),
  .errorhandling="pass") %dopar% {
    theta_mat = NULL
    var_mat = NULL
    CR_mat = NULL
    alpha_vec = NULL
    for(cnt in 1:length(rnames)){
      # RM = modelcases[cnt, 1]
      # OM = modelcases[cnt, 2]
      # VM = modelcases[cnt, 3]
      
      # tmptime1 <- Sys.time()
      
      n = nrow(nhanes0_comp) # n = 200
      # x = model.matrix(~ ., data = nhanes0_comp %>% select(!c(PEASCCT1, BPXPLS, BPXSY1, BPXSY2, BPXSY3, BPXDI1, BPXDI2, BPXDI3)))[, -1]
      x = model.matrix(~ ., data = nhanes0_comp %>% select(!c(BPXSY1, BPXSY2, BPXSY3, BPXDI1, BPXDI2, BPXDI3)))[, -1]
      # x = model.matrix(~ ., data = nhanes0_comp %>% select(!c(BPXSY1)))[, -1]
      x = scale(x)
      pcol = ncol(x)
      z = x
      # if(RM  == T){
      #   pi = 1 / (1 + exp(-(-.5 - x[,3] + 0.5 * x[,4] + 0.5 * x[,5] - 0.25 * x[,6])))
      # }else if(RM == F){
      #   #   pi = 1 / (1 + exp(-(- 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.75 * (x[,2] - 2)^2 -
      #   #                         0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
      #   # pi = ifelse(pi < 0.05, 0.05, pi)
      #   # pi = ifelse(pi > 0.7, 0.7, pi)
      #     # pi = 1 / (1 + exp(-(1 - 0.5 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
      #     #                       0.5 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2 )))
      #     # pi = pt(- 0.1 * (x[,1] - 2) * (x[,1] - 5) - 0.5 * (x[,2] - 2)^2 -
      #     #           0.2 * (x[,3] - 1.5)^2 * (x[,4] - 2)^2, 3)
      #   pi = 1 / (1 + exp(-(-.5 - 0.5 * (x[,3] - 2) * (x[,4] - 5) - 0.5 * (x[,4] - 2.5)^2 -
      #                         0.5 * (x[,5] - 1.5)^2 * (x[,6] - 2)^2 )))
      #     pi = ifelse(pi < 0.05, 0.05, pi)
      #     # pi = ifelse(pi > 0.7, 0.7, pi)
      #     summary(pi)
      #     mean(pi)
      # }

      vx = rep(1, n) # original
      
      y <- nhanes0_comp$BPXSY1
      theta = mean(nhanes0_comp$BPXSY1)
      
      delta = rbinom(n, 1, pi)
      Index_S = (delta == 1)
      
      theta_res = NULL
      var_res = NULL
      CR_res = NULL
      
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

      # ps.cv.RCAL <- glm.regu.cv(fold=3, nrho=1+20, y=delta, x=xR, loss="cal")
      # pihat.RCAL <- ps.cv.RCAL$sel.fit[,1]
      
      # ps.cv.RCAL <- glm.regu.cv(fold=5, nrho=1+10, y=delta, x=xR, loss="cal")
      # ps.cv.RCAL$sel.fit[,1] # pihat.cv.RCAL
      # ps.RCAL = glm.regu(delta, x1, rhos = rep(ps.cv.RCAL$sel.rho[1], ncol(x1)), loss = "cal") # Non default pihat.RCAL
      # ps.RCAL = glm.regu(delta, xR, rhos = rep(ps.cv.RCAL$sel.rho[2], ncol(xR)), loss = "cal") # ps.cv.RCAL$sel.rho[2] is default
      # ps.RCAL = glm.regu(delta, xR, rhos = rep(0.1, ncol(xR)), loss = "cal")
            
      # phihat.RCAL = c(ps.RCAL$inter, ps.RCAL$bet)
      
      # pihat.RCAL <- ps.RCAL$fit # pihat.cral
      # pihat.RCAL <- 1 / findphi0(phihat.RCAL, xR = xR, Z = xR, delta = delta, returnw = T)

      # sum(y_S / pihat.RCAL[Index_S]) / sum(1 / pihat.RCAL[Index_S]) #IPW.RCAL
      
      ps.cv.MLE <- glm.regu.cv(fold=3, nrho=1+20, y=delta, x=xR, loss="ml")
      pihat.MLE <- ps.cv.MLE$sel.fit[,1] # pihat using MLE
      # ps.MLE = glm.regu(delta, xR, rhos = rep(ps.cv.MLE$sel.rho[1], ncol(xR)), loss = "ml")
      
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

      or.cv <- glm.regu.cv(fold=3, nrho=1+20, y = y_S, x = xO_S, loss="gaus")
      yhat <- c( cbind(1,x)%*%or.cv$sel.bet[,1] )
      # if(OM){
      #   sel.rho = 0.6
      # }else if(OM){
      #   sel.rho = 0.2
      # } 
      # # or <- glm.regu(y = y_S, x = xO_S, rhos = rep(or.cv$sel.rho[1], ncol(xR)), loss="gaus")
      # or <- glm.regu(y = y_S, x = cbind(1, xO_S), rhos = rep(0.2, ncol(xR)+1), loss="gaus")
      # yhat <- c( cbind(1,xO) %*% or$bet )
      
      
      
      
      for(pimethod in c(1)){
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
        
        mn.aipw.res2 = mn.aipw(y, delta, fp=pihat, fo=rep(y_IPW, length(y)))
        
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
        var_res = c(var_res, IPW = mn.aipw.res2$var)  

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
        # for(type in c("EL", "GE")){
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
          orth_vec = -xR * ((1 - pihat) / pihat / gprime1(1 / pihat, type = type, del = del) * vx); 
          orth_vec = scale(orth_vec); orth_vec_S = orth_vec[Index_S,]

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
          res <- CVXR::solve(prob, solver = "MOSEK")
          if(res$status == "solver_error"){
            w_S = rep(NA, length(d_S))
          }else{
            w_S = drop(res$getValue(w))
          }
          type = ifelse(type == "GE", "HD", type)
          theta_res = c(theta_res, setNames(sum(y_S * w_S) / n, paste(type, 1, sep = "_"))) # EL1
          var_res = c(var_res, setNames(NA, paste(type, 1, sep = "_")))
          # var_res = c(var_res, setNames(sum((c(y_S * w_S + yhat[Index_S] * (1 - w_S),
          #                                      yhat[!Index_S]) - sum(y_S * w_S) / n)^2) / n / (n-1), paste(type, 0, sep = "")))

          tau1 = max(abs(t(xR_S) %*% w_S / n - colMeans(xR)));
          tau2 = max(abs(t(orth_vec_S) %*% w_S / n - colMeans(orth_vec)));
          
          # rho = .7
          for(rho in rho_vec){
            print(rho)
            if(type == "EL"){
              Phi_R <- Maximize(sum(log(w)))
            }else if(type == "ET"){
              Phi_R <- Maximize(sum(entr(w) - w))
            }else if(type == "GE"){
              Phi_R <- Maximize(sum(-power(w, del + 1) * del^(-1) * (del + 1)^(-1)))
            }
            
            constraints <- list(Z_St %*% w / n == Zbar,
                                abs(t(xR_S) %*% w / n - colMeans(xR)) <= rho * tau1,
                                abs(t(orth_vec_S) %*% w / n - colMeans(orth_vec)) <= rho * tau2)
            
            prob <- CVXR::Problem(Phi_R, constraints)
            
            
            res <- CVXR::solve(prob, solver = "MOSEK")
            # res <- CVXR::solve(prob, solver = "ECOS_BB")
            # installed_solvers()
            
            if(res$status == "solver_error"){
              w_S = rep(NA, length(d_S))
            }else{
              w_S = drop(res$getValue(w))
            }
            
            # max(abs(t(orth_vec_S) / pihat[Index_S] / n - colMeans(orth_vec)))
            
            # max(abs(t(xR_S) / pihat[Index_S] / n - colMeans(xR)))
            
            # max(abs(Z_St %*% w_S / n - Zbar));max(abs(t(xR_S) %*% w_S / n - colMeans(xR))); max(abs(t(orth_vec_S) %*% w_S / n - colMeans(orth_vec)));
            type = ifelse(type == "GE", "HD", type)
            theta_res = c(theta_res, setNames(sum(y_S * w_S) / n, paste(type, rho, sep = "_"))) # EL1
            
            var_res = c(var_res, setNames(sum((c(y_S * w_S + yhat[Index_S] * (1 - w_S),
                                                 yhat[!Index_S]) - sum(y_S * w_S) / n)^2) / n / (n-1), paste(type, rho, sep = "_")))
          }

        }

        CR_res = ifelse(abs(theta_res - theta) < qnorm(0.975) * sqrt(var_res), 1, 0)
      }
      
      # tmptime2 <- Sys.time()
      # tmptime2 - tmptime1
      
      nlmres = NULL
      alpha_vec = c(alpha_vec, nlmres$par)
      theta_mat = cbind(theta_mat, theta_res - theta)
      var_mat = cbind(var_mat, var_res)
      CR_mat = cbind(CR_mat, CR_res)
    }
    list(theta_mat, alpha = alpha_vec, var_mat, CR_mat)
  }
# final_res
# save.image(paste("rdata/", timenow0, ".RData", sep = ""))
# save.image(paste("res", ".RData", sep = ""))
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

boxres = t(do.call(cbind, lapply(final_res1, c)))[, -c(1,2)] + theta
colnames(boxres) <- row.names(final_res1[[1]])[-c(1,2)]
rho_vec1 = c(1, rho_vec)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Given:
# boxres  # matrix or data.frame with columns like EL_1, EL_0.95, ..., HD_0.35
# rho_vec1 <- c(1.00, 0.95, 0.75, 0.55, 0.35)

rho_levels <- as.numeric(rho_vec1)

# Make sure we have a data.frame and a row id for melting:
df <- as.data.frame(boxres) |>
  mutate(.row = row_number())

# Long format: split names into type (EL/ET/HD) and rho (numeric)
long <- df |>
  pivot_longer(-.row, names_to = "key", values_to = "value") |>
  separate(key, into = c("type","rho"), sep = "_", remove = TRUE) |>
  mutate(
    rho  = as.numeric(rho),
    rho  = factor(rho, levels = rho_levels, ordered = TRUE),
    type = factor(type, levels = c("EL","ET","HD"))
  )

plot_one <- function(tp, theta) {
  dat <- dplyr::filter(long, type == tp)
  
  # per-rho NA counts and a y-position just above the tallest observation
  lab_df <- dat %>%
    group_by(rho) %>%
    summarise(
      na_count = sum(is.na(value)),
      ymax     = suppressWarnings(max(value, na.rm = TRUE)),
      .groups  = "drop"
    )
  
  # vertical offset so labels sit above the boxes; handle all-NA edge case
  rng <- range(dat$value, na.rm = TRUE)
  off <- if (is.finite(diff(rng))) 0.03 * diff(rng) else 0.5
  lab_df$y <- ifelse(is.finite(lab_df$ymax), lab_df$ymax + off, theta + off)
  
  ggplot(dat, aes(x = rho, y = value)) +
    geom_boxplot(na.rm = TRUE) +
    geom_hline(yintercept = theta, color = "red", linewidth = 1) +
    geom_text(
      data = lab_df,
      aes(x = rho, y = y, label = na_count),
      color = "blue", size = 3.5
    ) +
    labs(title = paste("GEC estimation for", tp), x = expression(tau), y = "BPXSY1") +
    theme_minimal(base_size = 12) +
    coord_cartesian(clip = "off")  # show labels even if slightly above panel
}

p_EL <- plot_one("EL", theta)
p_ET <- plot_one("ET", theta)
p_HD <- plot_one("HD", theta)

# Print them (or save with ggsave)
p_EL; p_ET; p_HD


# for(ncols in 1:Ncol){
#   # boxres = t(do.call(cbind, lapply(final_res1, c)))[, ((ncols - 1) * Nrow + 1) : (ncols * Nrow)] + theta
#   boxres = t(do.call(cbind, lapply(final_res1, c)))[, -c(1,2)] + theta
#   colnames(boxres) <- row.names(final_res1[[1]])[-c(1,2)]
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

# final_res_var = lapply(final_res, function(x) x[[3]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
# res_var <- matrix(rowMeans(do.call(cbind, lapply(final_res_var, c)), na.rm = TRUE), 
#                ncol = ncol(final_res_var[[3]]))
# colnames(res_var) = rnames; row.names(res_var) = row.names(final_res_var[[3]])
# 
# na_res = matrix(colSums(apply(do.call(cbind, lapply(final_res_var, c)), 1, is.na)),
#                 ncol = ncol(final_res_var[[1]]))
# colnames(na_res) = rnames; row.names(na_res) = row.names(final_res1[[1]])
# na_res
# 
# (res_var - res3^2) / res3^2
# 
# final_res_CR = lapply(final_res, function(x) x[[4]])[sapply(final_res0, function(x) is.numeric(unlist(x)))]
# res_CR <- matrix(rowMeans(do.call(cbind, lapply(final_res_CR, c)), na.rm = TRUE), 
#                   ncol = ncol(final_res_CR[[4]]))
# colnames(res_CR) = rnames; row.names(res_CR) = row.names(final_res_CR[[4]])
# # res_CR
# 
# gsub("\\\\addlinespace", "", kable(cbind((res_var - res3^2) / res3^2, res_CR), "latex", booktabs = TRUE, digits = 3) %>%
#        kable_styling())

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Given:
# boxres  # matrix or data.frame with columns like EL_1, EL_0.95, ..., HD_0.35
# rho_vec1 <- c(1.00, 0.95, 0.75, 0.55, 0.35)

rho_levels <- as.numeric(rho_vec1)

# Make sure we have a data.frame and a row id for melting:
df <- as.data.frame(boxres) |>
  mutate(.row = row_number())

# Long format: split names into type (EL/ET/HD) and rho (numeric)
long <- df |>
  pivot_longer(-.row, names_to = "key", values_to = "value") |>
  separate(key, into = c("type","rho"), sep = "_", remove = TRUE) |>
  mutate(
    rho  = as.numeric(rho),
    rho  = factor(rho, levels = rho_levels, ordered = TRUE),
    type = factor(type, levels = c("EL","ET","HD"))
  )

plot_one <- function(tp, theta) {
  ggplot(filter(long, type == tp), aes(x = rho, y = value)) +
    geom_boxplot(na.rm = TRUE) +
    geom_hline(yintercept = theta, color = "red", linewidth = 1) +
    labs(title = tp, x = expression(rho), y = "boxres") +
    theme_minimal(base_size = 12)
}

p_EL <- plot_one("EL", theta)
p_ET <- plot_one("ET", theta)
p_HD <- plot_one("HD", theta)

# Print them (or save with ggsave)
p_EL; p_ET; p_HD

plot_rmse <- function(tp, theta) {
  
  dat <- dplyr::filter(long, type == tp)
  dat <- dplyr::filter(dat, rho < 0.55)
  # dplyr::filter(dat)
  
  # per-rho RMSE and NA counts
  rmse_df <- dat %>%
    group_by(rho) %>%
    summarise(
      rmse     = sqrt(mean((value - theta)^2, na.rm = TRUE)),
      na_count = sum(is.na(value)),
      .groups  = "drop"
    )
  
  # small vertical offset for labels
  y_max <- suppressWarnings(max(rmse_df$rmse, na.rm = TRUE))
  off   <- if (is.finite(y_max)) 0.03 * y_max else 0.5
  
  ggplot(rmse_df, aes(x = rho, y = rmse, group = 1)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    geom_text(
      aes(y = rmse + off, label = na_count),
      color = "blue", size = 3.5, na.rm = TRUE
    ) +
    labs(
      title = paste("GEC RMSE for", tp),
      x = expression(tau^"*"),
      y = "RMSE"
    ) +
    theme_minimal(base_size = 12) +
    coord_cartesian(clip = "off")
}
library(dplyr)
p_EL <- plot_rmse("EL", theta)
p_ET <- plot_rmse("ET", theta)
p_HD <- plot_rmse("HD", theta)

# Print them (or save with ggsave)
p_EL; p_ET; p_HD
