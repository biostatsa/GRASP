#' Perform simulation to compare different sample size determination methods
#'
#' This function generates pseudo-datasets to assess the performance of different
#' sample size determination (SSD) methods for estimating the required sample size for
#' prediction model development in both Study 1 and Study 2. It can compare the
#' effectiveness of models such as logistic regression (LR), linear discriminant
#' analysis (LDA), naive Bayes (NB), and support vector machine (SVM) with a linear
#' kernel, against the rules-of-thumb and Riley's method.
#' @param N The sample size of the underlying whole target population. Should be a
#' large number (e.g., N=100,000).
#' @param d The dimension of the feature vector X.
#' @param cov_mat The covariance matrix of the features. If missing, feature vector X
#' will be generated independently and identically from a normal distribution with mean 0
#' and standard deviation 1. If specified, feature vector X will be generated by a multivariate
#' normal distribution with mean 0 and covariance matrix cov_mat.
#' @param beta0 The coefficient value for the intercept of the underlying model that determines
#' the probability of each patient developing the disease status.
#' @param coeff Vector of the coefficient values for the feature vector X of the underlying model that determines
#' the probability of each patient developing the disease status. Dimension of coeff should be equal to d.
#' @param X_index The indices of the features assumed to be included in the prediction model.
#' @param diff_AUC The difference of the AUC value in the hypothesis test (i.e., theta1-theta0).
#' Default to 0.05.
#' @param alpha The significance level used to compute the (1-alpha)100% CI in Study 1 and the
#' desired type I error control in Study 2. Default to 0.05.
#' @param beta The desired type II error control in Study 2. Default to 0.2, which results in
#' the level of statistical power to be 0.8.
#' @param L The threshold of the length of the (1-alpha)100% CI. Default to 0.1.
#' @param seed The random seed value.
#' @param nBoots The number of bootstraps. Default to 500.
#' @param nCores Number of cores to use in computing results. Set to 1 to not use parallel computing.
#' Defaults to 1.
#' @param model Type of the prediction model. Can take value as "LR" (logistic regression), "LDA" (linear
#' distriminant analysis), "naive_bayes", and "svmLinear" (SVM with linear kernel). Other machine learning
#' (ML) algorithms incorporated in R package 'caret' can also be specified.
#' @returns `data_n` A matrix summarizing the mean and sd of the estimated sample size for
#' the diseased (n1), non-diseased (n0), and the entire (n) patient group based on different SSD methods.
#' @returns `data_L` A matrix summarizing the mean adn sd of the CI length based on different SSD methods.
#' @returns `data_power` A matrix summarizing the power of hypothesis test based on different SSD methods.
#' @returns `data_var` A matrix summarizing the true positive rate (TPR) and false positive rate (FPR) in
#' identifying the significant and insignificant variables in the prediction model. Only returned if
#' logistic regression model is fitted.
#' @export
sampsize_simulation <- function(N,
                                d,
                                cov_mat,
                                beta0,
                                coeff,
                                X_index,
                                diff_AUC=0.05,
                                alpha=0.05,
                                beta=0.2,
                                L=0.1,
                                seed,
                                nBoots=500,
                                nCores=1,
                                model=c("LR","LDA","naive_bayes","svmLinear")){

  ################################################################################
  set.seed(seed)

  if(length(coeff)!=d){
    stop("Dimension of the coefficient vector coeff should equal d.")
  }

  ################################################################################
  ## generate population-level data ##############################################
  if(missing(cov_mat)){
    X <- matrix(rnorm(N*d, mean=0, sd=1), nrow=N)
  }else{
    X <- MASS::mvrnorm(n = N, mu = numeric(d), Sigma = cov_mat)
  }

  eta_vec <- cbind(1,X) %*% c(beta0,coeff)
  prob_vec <- exp(eta_vec)/(1+exp(eta_vec))

  y <- rbinom(N, size=1, prob = prob_vec)

  data <- as.data.frame(cbind(X,y))
  colnames(data) <- c(paste0("X",seq(1,d,1)),"Y")
  data$Y <- ifelse(data$Y==1, "Disease", "Non_disease")
  data$Y <- factor(data$Y, levels = c("Non_disease","Disease"))

  ################################################################################
  ## randomly sample a dataset as D_prior to calculate AUC and R2_CSadj ##########
  h <- length(which(data$Y=="Non_disease"))/length(which(data$Y=="Disease"))
  pevent <- length(which(data$Y=="Disease"))/N

  n1_prior <- 200
  n0_prior <- ceiling(h*n1_prior)
  index1_prior <- sample(which(data$Y=="Disease"), n1_prior, replace = FALSE)
  index0_prior <- sample(which(data$Y=="Non_disease"), n0_prior, replace = FALSE)
  data_prior <- data[c(index1_prior,index0_prior),]

  data_remain <- data[-c(index1_prior,index0_prior),]

  ################################################################################
  ## create matrices to store results ############################################
  ## Study 1
  if(model=="LR"){
    data_n <- data.frame("group"=c("n1","n0","n"),
                         "rules_of_thumb"=numeric(3),
                         "Riley"=numeric(3),
                         "AUC_Study1"=numeric(3),
                         "AUC_Study2"=numeric(3))
    data_L <- data.frame("rules_of_thumb"=0,
                         "Riley"=0,
                         "AUC_Study1"=0)
    data_power <- data.frame("rules_of_thumb"=0,
                             "Riley"=0,
                             "AUC_Study2"=0)
    data_var <- data.frame("group"=c("TPR","FPR"),
                           "rules_of_thumb"=numeric(2),
                           "Riley"=numeric(2),
                           "AUC_Study1"=numeric(2),
                           "AUC_Study2"=numeric(2))
  }else{
    data_n <- data.frame("group"=c("n1","n0","n"),
                         "rules_of_thumb"=numeric(3),
                         "AUC_Study1"=numeric(3),
                         "AUC_Study2"=numeric(3))
    data_L <- data.frame("rules_of_thumb"=0,
                         "AUC_Study1"=0)
    data_power <- data.frame("rules_of_thumb"=0,
                             "AUC_Study2"=0)
  }

  ################################################################################
  ################################################################################
  p <- length(X_index)
  X_index_sig <- X_index[which(coeff[X_index]!=0)]

  formula <- paste0("Y~",colnames(data)[X_index[1]])
  for(i in 2:p){
    formula <- paste0(formula, "+", colnames(data)[X_index[i]])
  }
  formula <- as.formula(formula)

  ## calculate AUC* ##############################################################
  if(model=="LR"){
    model_prior <- glm(formula, data=data_prior, family="binomial")
    pred_prob_prior <- predict(model_prior, newdata=data_prior, type="response")
  }else if(model=="LDA"){
    model_prior <- MASS::lda(formula, data=data_prior)
    pred_prob_prior <- predict(model_prior, newdata=data_prior)$posterior[,2]
  }else{
    control <- caret::trainControl(method="repeatedcv",
                                   number=5,
                                   repeats=1,
                                   savePredictions=TRUE,
                                   classProbs=TRUE,
                                   summaryFunction=twoClassSummary)
    model_prior <- caret::train(formula,
                                data=data_prior,
                                trControl=control,
                                method=model,
                                metric="ROC")
    pred_prob_prior <- predict(model_prior, newdata=data_prior, type="prob")[,2]
  }

  AUC <- pROC::roc(data_prior$Y, pred_prob_prior, quiet = TRUE)$auc[1]

  if(AUC<=0.5){
    stop("Guess of the underlying AUC value <= 0.5.")
  }

  ## calculate R2_CSadj ##########################################################
  if(model=="LR"){
    log_Lnull <- length(which(data_prior$Y=="Disease"))*log(length(which(data_prior$Y=="Disease"))/nrow(data_prior))+
      (nrow(data_prior)-length(which(data_prior$Y=="Disease")))*log(1-length(which(data_prior$Y=="Disease"))/nrow(data_prior))
    log_Lmodel <- as.numeric(logLik(model_prior))
    LR <- -2*(log_Lnull-log_Lmodel)
    R2_CSapp <- 1-exp(-LR/nrow(data_prior))
    S_VH <- 1-p/LR
    R2_CSadj <- S_VH*R2_CSapp
  }

  ################################################################################
  ################################################################################
  ## estimate sample size needed based on different methods

  ## rules of thumb
  n1_rule_of_thumb <- p*10
  n0_rule_of_thumb <- ceiling(h*(n1_rule_of_thumb))

  if(model=="LR"){
    ## Riley
    data_n_Riley <- pmsampsize::pmsampsize(type = "b", rsquared = R2_CSadj, parameters = p, prevalence = pevent)
    n1_Riley <- ceiling((data_n_Riley$sample_size)*pevent)
    n0_Riley <- data_n_Riley$sample_size-n1_Riley
  }

  ## AUC_based
  ## Study 1
  data_n_AUC_s1 <- GRASP_L(theta=AUC,
                                 alpha=alpha,
                                 L=L,
                                 h=h)
  n1_AUC_s1 <- data_n_AUC_s1[1,2]
  n0_AUC_s1 <- data_n_AUC_s1[2,2]

  ## Study 2
  data_n_AUC_s2 <- GRASP_test(theta0=ifelse(AUC-diff_AUC>0.5, AUC-diff_AUC, 0.5),
                                    theta1=AUC,
                                    alpha=alpha,
                                    beta=beta,
                                    h=h)
  n1_AUC_s2 <- data_n_AUC_s2[1,2]
  n0_AUC_s2 <- data_n_AUC_s2[2,2]

  ################################################################################
  ################################################################################
  ## randomly sample based on estimated sample size
  sim_AUC <- function(data_remain, n1, n0){
    n1_index <- sample(which(data_remain$Y=="Disease"), n1, replace = FALSE)
    n0_index <- sample(which(data_remain$Y=="Non_disease"), n0, replace = FALSE)

    data_sim <- rbind(data_remain[n1_index,],data_remain[n0_index,])

    if(model=="LR"){
      model_sim <- glm(formula, data=data_sim, family="binomial")
      pred_prob_sim <- predict(model_sim, newdata=data_sim, type="response")
    }else if(model=="LDA"){
      model_sim <- MASS::lda(formula, data=data_sim)
      pred_prob_sim <- predict(model_sim, newdata=data_sim)$posterior[,2]
    }else{
      model_sim <- caret::train(formula,
                                data=data_sim,
                                trControl=control,
                                method = model,
                                metric ="ROC")
      pred_prob_sim <- predict(model_sim, newdata=data_sim, type="prob")[,2]
    }

    AUC_sim <- pROC::roc(data_sim$Y, pred_prob_sim, quiet = TRUE)$auc[1]

    if(model=="LR"){
      sig_index <- X_index[which(summary(model_sim)$coefficients[2:nrow(summary(model_sim)$coefficients),4]<0.05)]
      TPR <- length(which(sig_index %in% X_index_sig))/length(X_index_sig)
      FPR <- length(which(!(sig_index %in% X_index_sig)))/(p-length(X_index_sig))
    }

    if(model=="LR"){
      return(list("AUC_sim"=AUC_sim,
                  "TPR"=TPR,
                  "FPR"=FPR))
    }else{
      return(list("AUC_sim"=AUC_sim))
    }

  }

  ################################################################################
  ## parallel computing
  if(nCores>1){

    RNGkind("L'Ecuyer-CMRG")

    sim_results_rules_of_thumb <- parallel::mclapply(1:nBoots, function(b) sim_AUC(data_remain, n1_rule_of_thumb, n0_rule_of_thumb), mc.cores=nCores)
    sim_results_AUC_s1 <- parallel::mclapply(1:nBoots, function(b) sim_AUC(data_remain, n1_AUC_s1, n0_AUC_s1), mc.cores=nCores)
    sim_results_AUC_s2 <- parallel::mclapply(1:nBoots, function(b) sim_AUC(data_remain, n1_AUC_s2, n0_AUC_s2), mc.cores=nCores)

    sim_results_rules_of_thumb_AUC <- do.call(c, lapply(1:nBoots, function(i) sim_results_rules_of_thumb[[i]]$AUC_sim))
    sim_results_AUC_s1_AUC <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s1[[i]]$AUC_sim))
    sim_results_AUC_s2_AUC <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s2[[i]]$AUC_sim))

    if(model=="LR"){
      sim_results_Riley <- parallel::mclapply(1:nBoots, function(b) sim_AUC(data_remain, n1_Riley, n0_Riley), mc.cores=nCores)
      sim_results_Riley_AUC <- do.call(c, lapply(1:nBoots, function(i) sim_results_Riley[[i]]$AUC_sim))

      sim_results_rules_of_thumb_TPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_rules_of_thumb[[i]]$TPR))
      sim_results_Riley_TPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_Riley[[i]]$TPR))
      sim_results_AUC_s1_TPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s1[[i]]$TPR))
      sim_results_AUC_s2_TPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s2[[i]]$TPR))

      sim_results_rules_of_thumb_FPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_rules_of_thumb[[i]]$FPR))
      sim_results_Riley_FPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_Riley[[i]]$FPR))
      sim_results_AUC_s1_FPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s1[[i]]$FPR))
      sim_results_AUC_s2_FPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s2[[i]]$FPR))
    }
  }else{
    sim_results_rules_of_thumb <- lapply(1:nBoots, function(b) sim_AUC(data_remain, n1_rule_of_thumb, n0_rule_of_thumb))
    sim_results_AUC_s1 <- lapply(1:nBoots, function(b) sim_AUC(data_remain, n1_AUC_s1, n0_AUC_s1))
    sim_results_AUC_s2 <- lapply(1:nBoots, function(b) sim_AUC(data_remain, n1_AUC_s2, n0_AUC_s2))

    sim_results_rules_of_thumb_AUC <- do.call(c, lapply(1:nBoots, function(i) sim_results_rules_of_thumb[[i]]$AUC_sim))
    sim_results_AUC_s1_AUC <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s1[[i]]$AUC_sim))
    sim_results_AUC_s2_AUC <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s2[[i]]$AUC_sim))

    if(model=="LR"){
      sim_results_Riley <- lapply(1:nBoots, function(b) sim_AUC(data_remain, n1_Riley, n0_Riley))
      sim_results_Riley_AUC <- do.call(c, lapply(1:nBoots, function(i) sim_results_Riley[[i]]$AUC_sim))

      sim_results_rules_of_thumb_TPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_rules_of_thumb[[i]]$TPR))
      sim_results_Riley_TPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_Riley[[i]]$TPR))
      sim_results_AUC_s1_TPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s1[[i]]$TPR))
      sim_results_AUC_s2_TPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s2[[i]]$TPR))

      sim_results_rules_of_thumb_FPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_rules_of_thumb[[i]]$FPR))
      sim_results_Riley_FPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_Riley[[i]]$FPR))
      sim_results_AUC_s1_FPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s1[[i]]$FPR))
      sim_results_AUC_s2_FPR <- do.call(c, lapply(1:nBoots, function(i) sim_results_AUC_s2[[i]]$FPR))
    }
  }

  ################################################################################
  ## save results ################################################################
  ## Study 1
  data_n$rules_of_thumb[1] <- n1_rule_of_thumb
  data_n$rules_of_thumb[2] <- n0_rule_of_thumb
  data_n$rules_of_thumb[3] <- n0_rule_of_thumb+n1_rule_of_thumb
  data_n$AUC_Study1[1] <- n1_AUC_s1
  data_n$AUC_Study1[2] <- n0_AUC_s1
  data_n$AUC_Study1[3] <- n0_AUC_s1+n1_AUC_s1
  data_n$AUC_Study2[1] <- n1_AUC_s2
  data_n$AUC_Study2[2] <- n0_AUC_s2
  data_n$AUC_Study2[3] <- n0_AUC_s2+n1_AUC_s2

  data_L$rules_of_thumb <- quantile(sim_results_rules_of_thumb_AUC, probs = c(0.025,0.975))[2]-quantile(sim_results_rules_of_thumb_AUC, probs = c(0.025,0.975))[1]
  data_L$AUC_Study1 <- quantile(sim_results_AUC_s1_AUC, probs = c(0.025,0.975))[2]-quantile(sim_results_AUC_s1_AUC, probs = c(0.025,0.975))[1]

  theta0 <- ifelse(AUC-diff_AUC>0.5, AUC-diff_AUC, 0.5)
  Z_alpha <- qnorm(1-alpha/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

  data_power$rules_of_thumb <- ifelse((mean(sim_results_rules_of_thumb_AUC)-theta0)/sd(sim_results_rules_of_thumb_AUC)> Z_alpha, 1, 0)
  data_power$AUC_Study2 <- ifelse((mean(sim_results_AUC_s2_AUC)-theta0)/sd(sim_results_AUC_s2_AUC) > Z_alpha, 1, 0)

  if(model=="LR"){
    data_n$Riley[1] <- n1_Riley
    data_n$Riley[2] <- n0_Riley
    data_n$Riley[3] <- n0_Riley+n1_Riley

    data_L$Riley <- quantile(sim_results_Riley_AUC, probs = c(0.025,0.975))[2]-quantile(sim_results_Riley_AUC, probs = c(0.025,0.975))[1]

    data_power$Riley <- ifelse((mean(sim_results_Riley_AUC)-theta0)/sd(sim_results_Riley_AUC)> Z_alpha, 1, 0)

    data_var$rules_of_thumb[1] <- mean(sim_results_rules_of_thumb_TPR)
    data_var$rules_of_thumb[2] <- mean(sim_results_rules_of_thumb_FPR)
    data_var$Riley[1] <- mean(sim_results_Riley_TPR)
    data_var$Riley[2] <- mean(sim_results_Riley_FPR)
    data_var$AUC_Study1[1] <- mean(sim_results_AUC_s1_TPR)
    data_var$AUC_Study1[2] <- mean(sim_results_AUC_s1_FPR)
    data_var$AUC_Study2[1] <- mean(sim_results_AUC_s2_TPR)
    data_var$AUC_Study2[2] <- mean(sim_results_AUC_s2_FPR)
  }

  ################################################################################
  ################################################################################
  if(model=="LR"){
    return(list("data_n"=data_n,
                "data_L"=data_L,
                "data_power"=data_power,
                "data_var"=data_var))
  }else{
    return(list("data_n"=data_n,
                "data_L"=data_L,
                "data_power"=data_power))
  }

}




