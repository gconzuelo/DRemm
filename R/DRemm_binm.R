setwd(".")
packages <- c("data.table","tidyverse","haven","lubridate","here","VGAM","splines","parallel","foreach",
              "Formula","plotrix","TeachingDemos","plotmo","nnls","mvtnorm","SuperLearner","glmnet","tmle",
              "gam","earth","ranger","nnet","npcausal","rpart","tmle3","sl3")
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

cat("Check if library directory exists",'\n')
ifelse(!dir.exists(userLib), dir.create(userLib), FALSE)

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

dr_here()

##############################################################
#                        SuperLearner                        #
##############################################################
## SUPER LEARNER LIBRARY
## LOADING ALL
ranger_learner1  <- create.Learner("SL.ranger", tune=list(min.node.size=30, num.trees=500,max.depth=2))
ranger_learner2  <- create.Learner("SL.ranger", tune=list(min.node.size=30, num.trees=500,max.depth=3))
glmnet_learner   <- create.Learner("SL.glmnet", tune=list(alpha = seq(0,1,.2)))
earth_learner1   <- create.Learner("SL.earth",tune=list(degree=2))
earth_learner2   <- create.Learner("SL.earth",tune=list(degree=3))
gam_learner1     <- create.Learner("SL.gam",tune=list(deg.gam=2))
gam_learner2     <- create.Learner("SL.gam",tune=list(deg.gam=3))
gam_learner3     <- create.Learner("SL.gam",tune=list(deg.gam=4))
SL.library <- c("SL.glm", "SL.nnet", "SL.mean", ranger_learner1$names, ranger_learner1$names,
                glmnet_learner$names, earth_learner1$names, earth_learner2$names, gam_learner1$names,
                gam_learner2$names, gam_learner3$names)

## FINAL LIBRARY
sl.lib_mu <- sl.lib_pi <- SL.library

##############################################################
#                            SL 3                            #
##############################################################
ranger_lrn1 <- make_learner(Lrnr_ranger, min.node.size = 30, num.trees=500, max.depth=2)
ranger_lrn2 <- make_learner(Lrnr_ranger, min.node.size = 30, num.trees=500, max.depth=3)
glmnet_lrn  <- make_learner(Lrnr_glmnet, alpha = seq(0,1,.2))
earth_lrn1  <- make_learner(Lrnr_earth, degree = 2)
earth_lrn2  <- make_learner(Lrnr_earth, degree = 3)
gam_lrnr   <- Lrnr_pkg_SuperLearner$new("SL.gam")
glm_lrnr   <- Lrnr_pkg_SuperLearner$new("SL.glm")
mean_lrnr  <- Lrnr_pkg_SuperLearner$new("SL.mean")
nnet_lrnr  <- Lrnr_pkg_SuperLearner$new("SL.nnet")

# create a list with all learners
lrn_list <- list(ranger_lrn1, ranger_lrn2, glmnet_lrn, earth_lrn1, earth_lrn2, gam_lrnr, glm_lrnr, 
                 mean_lrnr, nnet_lrnr)

# define metalearners appropriate to data types
metalearner <- make_learner(Lrnr_nnls)

# Define the sl_Y and sl_A (we only need 1 because both are same type)
sl_lib <- Lrnr_sl$new(learners = lrn_list, metalearner = make_learner(Lrnr_nnls))

#####
## FUNCTION SET-UP
## TRUE VALUES
nsim  <- 1000
true1 <- 6
true2 <- 3

## DUMMY TABLE TO EXPORT THE RESULTS OF THE SIMULATION
cols  <- c("glm", "glm_strat","AIPW_strat", "TMLE_strat","AIPW_SS","TMLE_SS")
cols0 <- paste0(cols,"0"); cols1 <- paste0(cols,"1")
cols  <- c(cols0,cols1, "nsim")
res.est <- data.frame(matrix(nrow=nsim,ncol=length(cols)))
colnames(res.est) <- cols
res.se  <- res.est

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

## FUNCTION
interaction_sim <- function(counter, n, error_sd) {
  ##############################################################
  #                           SET-UP                           #
  ##############################################################
  i = counter
  n = n
  p = 2
  set.seed(i)
  
  ## CONFOUNDERS (C = 2)
  sigma <- matrix(0,nrow=p,ncol=p); diag(sigma) <- 1
  c     <- rmvnorm(n, mean=rep(0,p), sigma=sigma)
  
  # DESING MATRIX FOR THE OUTCOME MODEL
  muMatT <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  parms3 <- c(3,3)
  parms4 <- c(log(1.5),log(1.5))
  beta   <- parms3; beta <- c(120, beta)
  
  # DESIGN MATRIX FOR THE PROPENSITY SCORE MODEL
  piMatT  <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  piMatT2 <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  theta   <- c(-.5,parms4)
  theta2  <- c(-.5,parms4)
  mu      <- muMatT%*%beta
  
  # PROPENSITY SCORE MODEL
  pi   <- expit(piMatT%*%theta)
  pi_m <- expit(piMatT2%*%theta2)
  x    <- rbinom(n,1,pi)
  m    <- rbinom(n,1,pi_m)
  
  # OUTCOME MODEL: EXPOSURE VALUE UNDER M == 0 IS 6; VALUE UNDER M == 1 IS 3
  y <- x*6 + m*6 - 3*x*m + mu + rnorm(n,0,error_sd)
  
  # CORRECT MODEL SPECIFICATION
  # DATA
  dat <- data.frame(c,m,x,y); colnames(dat)[1:ncol(c)] <- paste("c",1:ncol(c),sep="")
  
  # ALL VARIABLES
  C <- muMatT[,-1];colnames(C) <- paste("C",1:ncol(muMatT[,-1]), sep="")
  X <- x
  M <- m
  Y <- y
  
  # C, X and Y STRATIFIED BY EFFECT MODIFIER
  C0 <- subset(C,M==0)
  C1 <- subset(C,M==1)
  X0 <- subset(X,M==0)
  X1 <- subset(X,M==1)
  Y0 <- subset(Y,M==0)
  Y1 <- subset(Y,M==1)
  
  ##############################################################
  #                    GLM & GLM stratified                    #
  ##############################################################
  
  # MODEL 1: GLM WITH INTERACTION
  mod1 <- lm(y~x+m+x*m+c1+c2, data=dat)
  
  # MODEL 2a, 2b: IN SUBPOPULATION M == 0 & M == 1
  mod2_0 <- lm(y~x+c1+c2, data=subset(dat,m==0))
  mod2_1 <- lm(y~x+c1+c2, data=subset(dat,m==1))
  
  # ESTIMATORS
  # POINT ESTIMATE
  # GLM
  res.est$glm0[i] <- coef(mod1)["x"]
  res.est$glm1[i] <- coef(mod1)["x"] + coef(mod1)["x:m"]
  
  # GLM STRATIFIED
  res.est$glm_strat0[i] <- coef(mod2_0)["x"]
  res.est$glm_strat1[i] <- coef(mod2_1)["x"]
  
  # STANDARD ERROR
  # CLOSED FORM OF THE SE
  # GLM
  res.se$glm0[i] <- sqrt(vcov(mod1)["x","x"])
  res.se$glm1[i] <- sqrt(vcov(mod1)["x","x"] + vcov(mod1)["m","m"] - 2*vcov(mod1)["x","m"])
  
  # GLM STRATIFIED
  res.se$glm_strat0[i] <- sqrt(vcov(mod2_0)["x","x"])
  res.se$glm_strat1[i] <- sqrt(vcov(mod2_1)["x","x"])
  
  ##############################################################
  #                       TMLE & AIPW                          #
  ##############################################################
  
  # TMLE for M == 0
  folds<-c(10)
  tmle_strat0 <- tmle(Y0,X0,C0, family="gaussian",
                      Q.SL.library=sl.lib_mu,g.SL.library=sl.lib_pi,
                      V=folds)
  
  # TMLE for M == 1
  tmle_strat1 <- tmle(Y1,X1,C1, family="gaussian",
                      Q.SL.library=sl.lib_mu,g.SL.library=sl.lib_pi,
                      V=folds)
  
  # Obtain PREDICTIONS from TMLE for M == 0 & 1 respectively
  pihat_0 <- tmle_strat0$g$g1W
  pihat_1 <- tmle_strat1$g$g1W
  
  pihat_0 <- ifelse(pihat_0 < quantile(pihat_0,0.025), quantile(pihat_0,0.025),pihat_0)
  pihat_0 <- ifelse(pihat_0 > quantile(pihat_0,0.975), quantile(pihat_0,0.975),pihat_0)
  
  pihat_1 <- ifelse(pihat_1 < quantile(pihat_1,0.025), quantile(pihat_1,0.025),pihat_1)
  pihat_1 <- ifelse(pihat_1 > quantile(pihat_1,0.975), quantile(pihat_1,0.975),pihat_1)
  
  # Obtain MUHAT for OVERALL and M == 0 (above) & 1 (under) respectively
  muhat_0   <- tmle_strat0$Qinit$Q[,2] * X0 + tmle_strat0$Qinit$Q[,1]*(1-X0)
  muhat_1   <- tmle_strat1$Qinit$Q[,2] * X1 + tmle_strat1$Qinit$Q[,1]*(1-X1)
  
  muhat1_0  <- tmle_strat0$Qinit$Q[,2]
  muhat1_1  <- tmle_strat1$Qinit$Q[,2]
  
  muhat0_0  <- tmle_strat0$Qinit$Q[,1]
  muhat0_1  <- tmle_strat1$Qinit$Q[,1]
  
  # ESTIMATORS
  # AIPW POINT ESTIMATE
  lower_bound <- min(y) - max(y)
  upper_bound <- max(y) - min(y)
  
  aipw_strat0 <- mean((((2*X0-1)*(Y0 - muhat_0))/((2*X0-1)*pihat_0 + (1-X0)) + muhat1_0 - muhat0_0))
  aipw_strat1 <- mean((((2*X1-1)*(Y1 - muhat_1))/((2*X1-1)*pihat_1 + (1-X1)) + muhat1_1 - muhat0_1))
  
  aipw_strat0 <- ifelse(aipw_strat0>upper_bound,upper_bound,aipw_strat0)
  aipw_strat0 <- ifelse(aipw_strat0<lower_bound,lower_bound,aipw_strat0)
  
  aipw_strat1 <- ifelse(aipw_strat1>upper_bound,upper_bound,aipw_strat1)
  aipw_strat1 <- ifelse(aipw_strat1<lower_bound,lower_bound,aipw_strat1)
  
  res.est$AIPW_strat0[i]  <- aipw_strat0
  res.est$AIPW_strat1[i]  <- aipw_strat1
  
  # TMLE POINT ESTIMATE
  res.est$TMLE_strat0[i]<- tmle_strat0$estimates$ATE$psi
  res.est$TMLE_strat1[i]<- tmle_strat1$estimates$ATE$psi
  
  # STANDARD ERROR
  # CLOSED FORM OF THE SE
  # TMLE
  res.se$TMLE_strat0[i] <- (1/(1-mean(M)))*sqrt(tmle_strat0$estimates$ATE$var.psi)
  res.se$TMLE_strat1[i] <- (1/(mean(M)))*sqrt(tmle_strat1$estimates$ATE$var.psi)
  
  # AIPW
  aipw_strat0_se <- sd((((2*X0-1)*(Y0 - muhat_0))/((2*X0-1)*pihat_0 + (1-X0)) + muhat1_0 - muhat0_0))/sqrt(nrow(C0))
  aipw_strat1_se <- sd((((2*X1-1)*(Y1 - muhat_1))/((2*X1-1)*pihat_1 + (1-X1)) + muhat1_1 - muhat0_1))/sqrt(nrow(C1))
  
  res.se$AIPW_strat0[i]  <- (1/(1-mean(M)))*aipw_strat0_se
  res.se$AIPW_strat1[i]  <- (1/(mean(M)))*aipw_strat1_se
  
  ##############################################################
  #                         AIPW w/SS                          #
  ##############################################################
  
  # AIPW with Sampling Splitting, Edward Kennedy npcausal
  aipw.SS0 <- ate(y = Y0, a = X0, x = C0, nsplits=10, sl.lib = sl.lib_pi)
  aipw.SS1 <- ate(y = Y1, a = X1, x = C1, nsplits=10, sl.lib = sl.lib_pi)
  res.est$AIPW_SS0[i] <- aipw.SS0$res$est[3]
  res.est$AIPW_SS1[i] <- aipw.SS1$res$est[3]
  res.se$AIPW_SS0[i] <- aipw.SS0$res$se[3]
  res.se$AIPW_SS1[i] <- aipw.SS1$res$se[3]
  
  ##############################################################
  #                         TMLE w/SS                          #
  ##############################################################
  
  # TMLE with Sample Splitting, TMLE3 (TLVERSE) package
  df0 <- as.data.frame(cbind(Y0, X0, C0)); names(df0) <- c("Y0", "X0", "C1", "C2")
  df1 <- as.data.frame(cbind(Y1, X1, C1)); names(df1) <- c("Y1", "X1", "C1", "C2")
  node_list0 <- list(W = c("C1", "C2"), A = "X0", Y = "Y0")
  node_list1 <- list(W = c("C1", "C2"), A = "X1", Y = "Y1")
  
  ate_spec <- tmle_ATE(treatment_level = 0, control_level = 1)
  learner_list <- list(A = sl_lib, Y = sl_lib)
  tmle_fit0 <- tmle3(ate_spec, df0, node_list0, learner_list)
  tmle_fit1 <- tmle3(ate_spec, df1, node_list1, learner_list)
  
  res.est$TMLE_SS0[i] <- abs(tmle_fit0$summary$tmle_est)
  res.est$TMLE_SS1[i] <- abs(tmle_fit1$summary$tmle_est)
  res.se$TMLE_SS0[i] <- tmle_fit0$summary$se
  res.se$TMLE_SS1[i] <- tmle_fit1$summary$se
  
  ##############################################################
  #                  PRINT AND UPDATE RESULTS                  #
  ##############################################################
  
  # COMPILING THE RESULTS
  res.est0 <- res.est %>% select(ends_with("0"))
  res.est1 <- res.est %>% select(ends_with("1"))
  
  res.se0  <- res.se %>% select(ends_with("0"))
  res.se1  <- res.se %>% select(ends_with("1"))
  
  res.cov0 <- res.est0-1.96*res.se0 < true1 & true1 < res.est0+1.96*res.se0
  res.cov1 <- res.est1-1.96*res.se1 < true2 & true2 < res.est1+1.96*res.se1
  
  res.width <- (res.est[-13]+1.96*res.se[-13]) - (res.est[-13]-1.96*res.se[-13])
  
  # FINAL DATASET TO BE EXPORTED
  tmp <- data.frame(rbind(c(n,apply(res.est0-true1,2,mean,na.rm=T), apply(res.est1-true2,2,mean,na.rm=T)),
                          c(n,apply((res.est0-true1)^2,2,mean,na.rm=T), apply((res.est1-true2)^2,2,mean,na.rm=T)),
                          c(n,apply(res.cov0,2,mean,na.rm=T), apply(res.cov1,2,mean,na.rm=T)),
                          c(n,apply(res.width,2,mean,na.rm=T)),
                          c(n,apply(res.se0,2,mean,na.rm=T),apply(res.se1,2,mean,na.rm=T))
                          ))
  tmp$nsim <- i
  
  rownames(tmp)      <- c("bias","rmse","cov","width","SE");
  colnames(tmp)[1]   <- "N"
  setDT(tmp, keep.rownames = TRUE)[]

  return(tmp)
}


## RUNNING THE FUNCTION
results <- lapply(1:nsim, function(x) interaction_sim(x,n=100,error_sd=3))

## COMPILING THE RESULTS
res <- do.call(rbind, results)
res <- as.data.frame(res)

write.table(res, file = "Results.txt", sep="\t", row.names = F, col.names = T)


