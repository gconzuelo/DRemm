setwd(".")
packages <- c("data.table","tidyverse","MASS","haven","lubridate","here","VGAM","splines","parallel","foreach",
              "Formula","plotrix","e1071","TeachingDemos","plotmo","nnls","mvtnorm","SuperLearner","glmnet","tmle",
              "gam","earth","ranger","nnet","npcausal","rpart","tmle3","sl3","mfp","boot","sandwich")
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

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

##############################################################
#                        SuperLearner                        #
##############################################################
ranger_learner  <- create.Learner("SL.ranger", tune=list(min.node.size = c(15, 30),
                                                         num.trees=c(1000),
                                                         max.depth=c(3,6)))
glmnet_learner  <- create.Learner("SL.glmnet", tune=list(alpha = seq(0,1,.2)))
svm_learner     <- create.Learner("SL.svm",tune=list(nu = c(.5,.75),degree=3))
earth_learner   <- create.Learner("SL.earth",tune=list(degree=c(2,3)))
gam_learner     <- create.Learner("SL.gam",tune=list(deg.gam=c(2,3,4)))
SL.library <- c("SL.glm","SL.mean","SL.glm.interaction","SL.bayesglm", ranger_learner$names,
                glmnet_learner$names, svm_learner$names, earth_learner$names, gam_learner$names)

## FINAL LIBRARY
sl.lib_mu <- sl.lib_pi <- SL.library

##############################################################
#                           FUNCTION                         #
##############################################################

nsim  <- 200
emm.cont <- function(counter, func_type) {
  ##############################################################
  #                           SET-UP                           #
  ##############################################################
  i = counter
  n = 500
  p = 2
  mu.p <- 0
  sd.p <- 1
  set.seed(888 + i)
  
  ## Covariance Matrix
  sigma <- matrix(0, nrow=p, ncol=p)
  diag(sigma) <- sd.p^2
  #
  c <- rmvnorm(n, mean = rep(mu.p, p), sigma = sigma)
  
  # Outcome Model Matrix
  muMat <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  parms <- c(3,3)
  beta  <- c(120, parms)
  
  # Treatment Model Matrix (piMat and muMat are not necessarily the same always)
  piMat <- muMat
  parmT <- c(log(1.5),log(1.5))
  theta <- c(-.5, parmT)
  
  # Effect Modifier Model Matrix (we restricted to positive values only)
  theta2 <- c(15, parms)
  mu_m <- pmax(0, piMat%*%theta2)
  
  # For each individual: mu = mean(Y); pi = Pr(trt)
  # Notice: expit(-.5) ~ mean(pi)
  mu <- muMat%*%beta
  pi <- expit(piMat%*%theta)
  x  <- rbinom(n,1,pi)
  m  <- (runif(n, min = 0, max = max(mu_m))*13)/max(mu_m)
  
  # Outcome Model
  yq <- y <- x*6 + m*2.5 + x*(m-6)^2 + mu + rnorm(n,0,6)
  yl <- x*6 + m*2.5 + x*log(m) + mu + rnorm(n,0,6)
  yc <- x*6 + m*2.5 + x*(4*sqrt(9*m)*as.numeric(m<2) + as.numeric(m>=2)*(abs(m-6)^(2))) + mu + rnorm(n,0,6)
  
  if (func_type == "quadratic") {
    y <- yq
  } else if (func_type == "incmono") {
    y <- yl
  } else if (func_type == "complex"){
    y <- yc
  }

  # Data
  dat <- data.frame(c,m,x,y); colnames(dat)[1:ncol(c)] <- paste("c",1:ncol(c),sep="")
  
  C <- muMat[,-1]; colnames(C) <- paste("C",1:ncol(muMat[,-1]), sep="")
  X <- x
  M <- m
  Z <- cbind(M,C)
  Y <- y
  
  m.value <- seq(0.1, 13.0, by = .1)
  
  folds<-c(10)
  tmle_res <- tmle(Y, X, Z, family = "gaussian",
                   Q.SL.library = sl.lib_mu, g.SL.library = sl.lib_pi,
                   V = folds)

  # PiHat from TMLE, including boundaries
  pihat <- tmle_res$g$g1W
  pihat <- ifelse(pihat < quantile(pihat,0.025), quantile(pihat,0.025),pihat)
  pihat <- ifelse(pihat > quantile(pihat,0.975), quantile(pihat,0.975),pihat)

  # Muhat from TMLE
  muhat      <- tmle_res$Qinit$Q[,2] * X + tmle_res$Qinit$Q[,1] * (1-X)
  muhat_star <- tmle_res$Qstar[,2] * X + tmle_res$Qstar[,1] * (1-X)
  muhat1     <- tmle_res$Qinit$Q[,2]
  muhat0     <- tmle_res$Qinit$Q[,1]

  ##############################################################
  #                       COMPUTE AIPW                         #
  ##############################################################
  # Obtaining the EFF
  aipw_EFF <- as.numeric((((2*X-1) * (Y - muhat)) / ((2*X-1) * pihat + (1-X)) + muhat1 - muhat0))
  X = cbind(0, M)

  # Obtaining the estimates
  aipw_eff <- SuperLearner(Y = aipw_EFF,
                           X = as.data.frame(X),
                           family = gaussian,
                           SL.library=sl.lib_mu,
                           method = "method.NNLS")

  mu_aipw_eff <- predict(aipw_eff, newdata = data.frame(V1=0, M=m.value), onlySL = T)$pred

  ##############################################################
  #                      COMPUTE SPLINE                        #
  ##############################################################
  
  ## Spline model
  mspl <- glm(y ~ x + bs(m, df = 3, degree = 3) + bs(x*m, df = 3, degree = 3) + c1 + c2,
              data = dat, family = gaussian(link = "identity"))
  
  # Predictions and difference calculation
  y_smooth0  <- predict(mspl, newdata = data.frame(x = rep(0,130), m = m.value, c1 = rep(0,130), c2 = rep(0,130)), type = "response")
  y_smooth1  <- predict(mspl, newdata = data.frame(x = rep(1,130), m = m.value, c1 = rep(0,130), c2 = rep(0,130)), type = "response")
  #
  mu_spline  <- y_smooth1 - y_smooth0
  
  ##############################################################
  #                       COMPUTE MFP                          #
  ##############################################################
  
  ## Frac. polynomials model
  fp <- mfp(y ~ x + fp(m, df = 4) + fp(x*m, df = 4) + c1 + c2, data = dat, family = gaussian(link = "identity"))
  
  # Predictions and difference calculation
  y_fp0  <- predict(fp, newdata = data.frame(x = rep(0,130), m = m.value, c1 = rep(0,130), c2 = rep(0,130)), type = "response")
  y_fp1  <- predict(fp, newdata = data.frame(x = rep(1,130), m = m.value, c1 = rep(0,130), c2 = rep(0,130)), type = "response")
  mu_fp  <- y_fp1 - y_fp0
  
  tmp <- data.frame(sim = i, psi_aipw = mu_aipw_eff, psi_spline = mu_spline, psi_fp = mu_fp, m = m.value)
  
  return(tmp)
  
}

## RUNNING THE FUNCTION
## func_type: quadratic == "quadratic"; increasing monotonic == "incmono"; complex == "complex
res <- lapply(1:nsim, function(x) {emm.cont(x, func_type = "complex")})

## RESULTS
results <- do.call(rbind, res)

write.table(results, file = "Results.txt", sep="\t", row.names = F, col.names = T)
