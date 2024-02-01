#' Estimate the required sample size based on AUC for Study 1
#'
#' This function calculates the lower and upper bounds of the required sample size
#' for the prediction model development to ensure that the length of the (1-alpha)100%
#' confidence interval (CI) for the empirical AUC estimator remains within a specific threshold.
#' @param theta The guess of the AUC value that the prediction model can potential achieve.
#' @param h The ratio of non-diseased to diseased patients, and we always assume that h >= 1.
#' @param alpha The significance level used to compute the (1-alpha)100% CI. Default to 0.05,
#' which results in a 95% CI.
#' @param L The threshold of the length of the (1-alpha)100% CI. Default to 0.1.
#' @returns `sampsize` A matrix summarizing the lower and upper bounds of the required sample size
#' for the diseased (n1), non-diseased (n0), and the entire (n) patient group.
#' @export
GRASP_L <- function(theta,
                          h,
                          alpha=0.05,
                          L=0.1){

  ################################################################################
  if(missing(theta)){
    stop("Error: guess of the underlying AUC value must be given.")
  }

  if(theta<=0 | theta>=1){
    stop("Error: guess of the underlying AUC value must be between 0 and 1.")
  }

  if(missing(h)){
    stop("Error: the ratio of non-disease to disease patients must be given.")
  }

  if(h<1){
    stop("Error: the ratio of non-disease to disease patients must be greater than or equal to 1.")
  }

  if(alpha<=0 | alpha>=1){
    stop("Error: the significance level must be between 0 and 1.")
  }

  if(L<=0 | L>=1){
    stop("Error: the threshold of the length of the CI must be between 0 and 1.")
  }



  ################################################################################
  ## B1, B2 functions ############################################################
  ## n0: non-disease patients
  ## n1: disease patients

  B1 <- function(n1,theta){
    n1-(n1-1)^2/(12*(h*n1-1)*theta*(1-theta))
  }

  B2 <- function(n1,theta){
    1-(h*n1+n1-2)*min(theta,1-theta)/max(theta,1-theta)+4/(3*max(theta,1-theta))*sqrt(2*min(theta,1-theta)*(h*n1-1)*(n1-1))
  }

  ################################################################################
  Z_alpha <- qnorm(1-alpha/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  V_cutoff <- L^2/(4*Z_alpha^2)

  ################################################################################
  ## calculate the upperbound of n1 ##############################################
  n1_u <- ceiling(theta*(1-theta)/V_cutoff)

  ################################################################################
  ## find the lowerbound of n1 ###################################################
  data_diff <- data.frame("n1"=seq(1,1000,1),
                          "diff"=numeric(length(seq(1,1000,1))))

  for(i in 1:nrow(data_diff)){
    if((data_diff$n1[i]-1)/(h*data_diff$n1[i]-1)<=2*min(theta,1-theta)){
      var_l <- theta*(1-theta)*B1(data_diff$n1[i],theta)/(h*data_diff$n1[i]^2)
    }else(
      var_l <- theta*(1-theta)*B2(data_diff$n1[i],theta)/(h*data_diff$n1[i]^2)
    )

    data_diff$diff[i] <- var_l-V_cutoff
  }

  if(length(which(data_diff$diff>0))>0 & length(which(data_diff$diff<=0))>0){
    n1_l <- data_diff$n1[which(data_diff$diff<=0)[which.max(data_diff$diff[which(data_diff$diff<=0)])]]
  }else{
    stop("Change point is not identified.")
  }

  ################################################################################
  ## return results ##############################################################
  sampsize <- data.frame("group"=c("disease:n1","non-disease:n0","all:n"),
                         "lowerbound"=c(n1_l, ceiling(h*n1_l), ceiling(h*n1_l)+n1_l),
                         "upperbound"=c(n1_u, ceiling(h*n1_u), ceiling(h*n1_u)+n1_u))

  return(sampsize)
}

#' Estimate the required sample size based on AUC for Study 2
#'
#' This function calculates the lower and upper bounds of the required sample size
#' for the prediction model development to test hypothesis with a desired type I error
#' control and a predetermined level of statistical power.
#' @param theta0 The AUC value specified in the null hypothesis, H0: AUC=theta0.
#' @param theta1 The AUC value specified in the alternative hypothesis, Ha: AUC=theta1>thata0.
#' @param h The ratio of non-diseased to diseased patients, and we always assume that h >= 1.
#' @param alpha The desired type I error control. Default to 0.05.
#' @param beta The desired type II error control. Default to 0.2, which results in
#' the level of statistical power to be 0.8.
#' @returns `sampsize` A matrix summarizing the lower and upper bounds of the required sample size
#' for the diseased (n1), non-diseased (n0), and the entire (n) patient group.
#' @export
GRASP_test <- function(theta0,
                             theta1,
                             h,
                             alpha=0.05,
                             beta=0.2){

  ################################################################################
  if(missing(theta0)){
    stop("Error: the AUC value in the null hypothesis must be given.")
  }

  if(missing(theta1)){
    stop("Error: the AUC value in the alternative hypothesis must be given.")
  }

  if(theta0<=0 | theta0>=1){
    stop("Error: the AUC value in the null hypothesis must be between 0 and 1.")
  }

  if(theta1<=0 | theta1>=1){
    stop("Error: the AUC value in the alternative hypothesis must be between 0 and 1.")
  }

  if(theta1<=theta0){
    stop("Error: the AUC value in the alternative hypothesis must be greater than the AUC value in the null hypothesis.")
  }

  if(missing(h)){
    stop("Error: the ratio of non-disease to disease patients must be given and we always assume h>=1.")
  }

  if(h<1){
    stop("Error: the ratio of non-disease to disease patients must be greater than or equal to 1.")
  }

  if(alpha<=0 | alpha>=1){
    stop("Error: the desired type I error control must be between 0 and 1.")
  }

  if(beta<=0 | beta>=1){
    stop("Error: the desired type II error control must be between 0 and 1.")
  }

  ################################################################################
  ## B1, B2 functions ############################################################
  ## n0: non-disease patients
  ## n1: disease patients

  B1 <- function(n1,theta){
    n1-(n1-1)^2/(12*(h*n1-1)*theta*(1-theta))
  }

  B2 <- function(n1,theta){
    1-(h*n1+n1-2)*min(theta,1-theta)/max(theta,1-theta)+4/(3*max(theta,1-theta))*sqrt(2*min(theta,1-theta)*(h*n1-1)*(n1-1))
  }

  ################################################################################
  ## calculate the upperbound of n1 ##############################################
  Z_alpha <- qnorm(1-alpha/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  Z_beta <- qnorm(beta, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

  n1_u <- ceiling(((Z_beta*sqrt(theta1*(1-theta1))-Z_alpha*sqrt(theta0*(1-theta0)))/(theta0-theta1))^2)

  ################################################################################
  ## find the lowerbound of n1 ###################################################
  data_diff <- data.frame("n1"=seq(1,1000,1),
                          "diff"=numeric(length(seq(1,1000,1))))

  for(i in 1:nrow(data_diff)){
    if((data_diff$n1[i]-1)/(h*data_diff$n1[i]-1)<=2*min(theta0,1-theta0)){
      var0_l <- theta0*(1-theta0)*B1(data_diff$n1[i],theta0)/(h*data_diff$n1[i]^2)
    }else(
      var0_l <- theta0*(1-theta0)*B2(data_diff$n1[i],theta0)/(h*data_diff$n1[i]^2)
    )

    if((data_diff$n1[i]-1)/(h*data_diff$n1[i]-1)<=2*min(theta1,1-theta1)){
      var1_l <- theta1*(1-theta1)*B1(data_diff$n1[i],theta1)/(h*data_diff$n1[i]^2)
    }else(
      var1_l <- theta1*(1-theta1)*B2(data_diff$n1[i],theta1)/(h*data_diff$n1[i]^2)
    )

    data_diff$diff[i] <- (Z_alpha*sqrt(var0_l)+theta0-theta1)/sqrt(var1_l)-Z_beta

  }

  if(length(which(data_diff$diff>0))>0 & length(which(data_diff$diff<=0))>0){
    n1_l <- data_diff$n1[which(data_diff$diff<=0)[which.max(data_diff$diff[which(data_diff$diff<=0)])]]
  }else{
    stop("Change point is not identified.")
  }

  ################################################################################
  ## return results ##############################################################
  sampsize <- data.frame("group"=c("disease:n1","non-disease:n0","all:n"),
                         "lowerbound"=c(n1_l, ceiling(h*n1_l), ceiling(h*n1_l)+n1_l),
                         "upperbound"=c(n1_u, ceiling(h*n1_u), ceiling(h*n1_u)+n1_u))

  return(sampsize)

}



