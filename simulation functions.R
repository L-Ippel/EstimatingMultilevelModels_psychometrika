##############
## simulation 2 revision psychometrika
##############

#install.packages('mvtnorm')
library(mvtnorm)
#note if you wish to use unstandardized varables make sure that
# the number of variances given equals the number of variables

genData <- function(n              = 50000,
                    j              = 2000,
                    n.level2       = 6, #this includes intercept and dummy variabels
                    #(hard coded: 2 cat. variables 1 with 2 cat, 1 with 3 cat)
                    n.level1       = 8,
                    # (hard coded 1 cat. variable with 4 categories)
                    fixed.coef     = c(100, runif(13, 0, 3)),
                    random.var     = runif(5, 2, 10),
                    cor.random.var = .5,
                    resid.var      = 5,
                    collinearity   = .15,
                    mean.level2    = 0,
                    mean.level1    = 0,
                    sd.level2      = 1,
                    sd.level1      = 1){
  level2Data <- genLevel2Data(j,
                              n.level2,
                              collinearity,
                              meanX = mean.level2,
                              var   = sd.level2**2)
  level1Data <- genLevel1Data(n,
                              meanX = mean.level1,
                              sdX   = sd.level1,
                              n.level1,
                              collinearity)

  id <- sample(1:j, n, replace = TRUE)
  if(length(random.var) != 1){
    randomCovar <- genCovarMatrix(cor.random.var, varX = random.var)
    randomCoef  <- rmvnorm(j, sigma = randomCovar)
  }
  else{
    randomCoef <- as.matrix(rnorm(j, 0, sqrt(random.var)) )
  }
  coefMatrix <- cbind(randomCoef[id, ] + rep(fixed.coef[1:ncol(randomCoef)], each = n),
                      matrix(
                        rep(fixed.coef[-(1:ncol(randomCoef))], each = n),
                      nrow = n))
  dataSet   <- as.data.frame(cbind(id, 1, level1Data, level2Data[id,]))
  dataSet$y <- rowSums(dataSet[, -1] * coefMatrix) + rnorm(n, 0, sqrt(resid.var))
  dataSet   <- as.data.frame(cbind(id  = dataSet[, 1],  #id
                                   y   = dataSet$y,
                                   int = dataSet[, 2], #intercept
                                   level1Data))
  temp.dim <- dim(dataSet)[2]
  for(b in 4:(temp.dim - 3)){
    names(dataSet)[b] <- paste0('level1.x', b - 3, sep = '')
  }
  names(dataSet)[(temp.dim - 2):temp.dim] <- c('night',
                                               'morning',
                                               'afternoon')
  dataSet <- as.data.frame(cbind(dataSet, level2Data[id,]))

  for(b in (1 + temp.dim):(dim(dataSet)[2] - 3)){
    names(dataSet)[b] <- paste0('level2.x', b - temp.dim, sep = '')
  }

  names(dataSet)[-c(1:(dim(dataSet)[2] - 3))] <- c('sex', 'edu.low', 'edu.high')

  return(dataSet)
}

genCovarMatrix <-function(cor.random.var,
                         varX){
  covar.level2 <- matrix(NA, nrow = length(varX), ncol = length(varX))
  for(i in 1:(length(varX))){
    for(j in 1:length(varX)){
      if(i == j){
        covar.level2[i, j] <- varX[i]
      }
      else{
        covar.level2[i, j] <- covar.level2[j, i] <-
          cor.random.var * sqrt(varX[i]) * sqrt(varX[j])
      }
    }
  }
  return(covar.level2)
}

genLevel2Data <- function(j,
                          n.level2,
                          collinearity,
                          meanX = 0,
                          varX = 1){
  if(length(meanX) == 1){
    meanX <- rep(meanX, n.level2 - 4) # 4= intercept, dummy 'sex' + 2 dummies 'low, middle'
  }
  if(length(varX) != 1) {
    covarMat <- genCovarMatrix(cor.random.var = collinearity, var = varX)
  }
  else{
    covarMat      <- matrix(collinearity, nrow = n.level2-4,
                             ncol = n.level2 - 4)
    diag(covarMat) <- 1
  }
  catVar <- matrix(NA, j, 2)
  for(i in 1:j){
    catVar[i, ] <- sample(c(0, 0, 1), 2, replace = F)
  }
  DataMatrix <- cbind(rmvnorm(j, mean = meanX, sigma = covarMat),
    rbinom(j, 1, .5), # dummy 1 represents gender
    catVar) # categorical variable with 3 categories, 2 dummies

  return(DataMatrix)
}

genLevel1Data <- function(n,
                          meanX,
                          sdX,
                          n.level1,
                          collinearity){
  xMat <- matrix(NA, nrow = n, ncol = 3)
  for(i in 1:n){
    xMat[i,] <- sample(c(0, 0, 0, 1), 3, replace = F)
  }
  if(length(meanX) != n.level1){
    meanX <- rep(meanX, n.level1 - 3) #hardcoded level 1 cat. variable with 4 categoriesc
  }
  if(length(sdX) != n.level1){
    sdX <- rep(sdX, n.level1 - 3)
  }
  sigma <- genCovarMatrix(cor.random.var = collinearity, var = sdX**2)
  xMat  <- cbind(rmvnorm(n, meanX, sigma), xMat)
  return(xMat)
}

#offline EM algorithm
# emAlgorithm <- function(data         = genData(),
#                         J            = NULL,
#                         start.fixed  = c(90, rnorm(13)),
#                         start.random = runif(5, 0, 4),
#                         data.random  = c(3, 9:12),
#                         start.cor    = .15,
#                         start.res    = 2,
#                         max.iter     = 20,
#                         crit.value   = .0001){
#   if(is.null(J)) J <- length(unique(data$id))
#   if(!is.matrix(start.random)){
#     start.random <- genCovarMatrix(cor.random.var = start.cor,
#                                    varX           = start.random)
#   }
#   n        <- nrow(data)
#   theta_j  <- list()
#
#   Xmat.inv <- solve(t(as.matrix(data[, -c(1,2)])) %*% as.matrix(data[, -c(1, 2)]))
#   XY.vec   <- t(as.matrix(data[, -c(1, 2)]))%*% as.matrix(data[, 2])
#   id.statistics <- unique(data$id)
#
#   for(t in 1:J){
#     temp             <- list()
#     temp$id          <- id.statistics[t]
#     temp$n_j         <- nrow(data[data$id == id.statistics[t], ])
#     temp$z_sq        <- matrix(0, nrow = length(data.random), ncol = length(data.random)) + t(as.matrix(
#                           data[data$id == id.statistics[t], data.random]))%*%
#                           as.matrix(data[data$id == id.statistics[t], data.random])
#
#     temp$zx_mat      <- matrix(0, nrow = nrow(Xmat.inv), ncol = length(data.random)) +
#                           as.matrix(t(data[data$id == id.statistics[t], -c(1, 2)])) %*%
#                           as.matrix(data[data$id == id.statistics[t], data.random])
#     temp$y_sq        <- as.matrix(t(data$y[data$id == id.statistics[t]])) %*%
#                           as.matrix(data$y[data$id == id.statistics[t]])
#     temp$x_sq        <- as.matrix(t(
#                           data[data$id == id.statistics[t], -c(1, 2)])) %*%
#                           as.matrix(data[data$id == id.statistics[t], -c(1, 2)])
#     temp$xy          <- as.matrix(t(data$y[data$id == id.statistics[t]])) %*%
#                           as.matrix(data[data$id == id.statistics[t], -c(1,2)])
#
#     temp$zy          <- t(as.matrix(data[data$id == id.statistics[t], data.random])) %*%
#                           as.matrix(data[data$id == id.statistics[t], 2])
#
#
#     theta_j[[t]]     <- temp
#   }
#     for(i in 1:max.iter){
#       cat(i,"\n")
#       T1         <- rep(0, nrow(Xmat.inv))
#       T2         <- matrix(0, nrow = length(data.random), ncol = length(data.random))
#       T3         <- 0
#       random.inv <- solve(start.random)
#       for(s in 1:J){
#         theta_j[[s]]$c_inv <- solve( theta_j[[s]]$z_sq + start.res * random.inv)
#
#         theta_j[[s]]$random_coef <- theta_j[[s]]$c_inv %*% (theta_j[[s]]$zy - t(theta_j[[s]]$zx_mat)
#                                                                                    %*% as.matrix(start.fixed))
#         theta_j[[s]]$t1_j <- theta_j[[s]]$zx_mat %*% theta_j[[s]]$random_coef
#         theta_j[[s]]$t2_j <- theta_j[[s]]$random_coef %*% t(theta_j[[s]]$random_coef) +
#                                start.res * theta_j[[s]]$c_inv
#         theta_j[[s]]$t3_j <- t((data$y[data$id == id.statistics[s]] -
#                               as.matrix(data[data$id == id.statistics[s], -c(1, 2)]) %*%
#                               start.fixed -
#                                 as.matrix(data[data$id == id.statistics[s], data.random]) %*%
#                                 theta_j[[s]]$random_coef)) %*%
#                               (data$y[data$id == id.statistics[s]] - as.matrix(
#                               data[data$id == id.statistics[s], -c(1, 2)] ) %*%
#                               start.fixed -as.matrix(
#                               data[data$id == id.statistics[s], data.random]) %*%
#                               theta_j[[s]]$random_coef) + start.res * sum(diag(theta_j[[s]]$c_inv  %*%
#                               theta_j[[s]]$z_sq))
#         T1 <- T1 + theta_j[[s]]$t1_j
#         T2 <- T2 + theta_j[[s]]$t2_j
#         T3 <- T3 + theta_j[[s]]$t3_j
#       }
#       fixed  <- Xmat.inv %*% (XY.vec - T1)
#       random <- T2 / J
#       res    <- as.numeric(T3 / n)
#       if(abs(max(fixed - start.fixed))   < crit.value &
#          abs(max(random - start.random)) < crit.value &
#          abs(res - start.res)            < crit.value) break
#       else{
#         cat(abs(max(fixed - start.fixed)), '\n')
#         cat(abs(max(random - start.random)), '\n')
#         cat(abs(res - start.res), '\n')
#
#         start.fixed  <- fixed
#         start.random <- random
#         start.res    <- res
#       }
#     }
#   theta <- list()
#   theta$fixed_coef_hat  <- fixed
#   theta$t1	            <- T1
#
#   theta$random_var_hat	<- random
#   theta$t2		          <- T2
#
#   theta$resid_var_hat	  <- res
#   theta$t3              <- T3
#
#   theta$n	              <- n
#   theta$j	              <- J
#
#   theta$xy_vector	      <- XY.vec
#   theta$x_inv	          <- Xmat.inv
#
#
#   return(list(model   = theta,
#               unit    = theta_j,
#               id_list = id.statistics))
# }

try.record <- function(id, id.param) {
  tryCatch(
    id.param[[id]],
    error = function(e)
    {
      return(NULL)
    }
  )
}

semaUpdate <- function(theta_jList,
                       theta,
                       J,
                       n) {
  theta$t1                 <- rep(0, length(theta$t1))
  theta$t2                 <-
    matrix(0, nrow = nrow(theta$t2), ncol = ncol(theta$t2))
  theta$t3                 <- 0
  random.inv               <- solve(theta$random_var_hat)

  for (i in 1:J) {
    theta_jList[[i]]$c_inv       <- solve(theta_jList[[i]]$z_sq +
                                            as.numeric(theta$resid_var_hat) * random.inv)
    theta_jList[[i]]$random_coef <- theta_jList[[i]]$c_inv %*%
      (
        as.matrix(theta_jList[[i]]$zy) -
          t(theta_jList[[i]]$zx_mat) %*%
          as.matrix(theta$fixed_coef_hat)
      )

    theta_jList[[i]]$t1_j        <- theta_jList[[i]]$zx_mat %*%
      as.matrix(theta_jList[[i]]$random_coef)

    theta_jList[[i]]$t2_j        <- theta_jList[[i]]$random_coef %*%
      t(theta_jList[[i]]$random_coef) +
      as.numeric(theta$resid_var_hat) * theta_jList[[i]]$c_inv

    theta_jList[[i]]$t3_j        <- theta_jList[[i]]$y_sq +
      t(theta$fixed_coef_hat) %*% theta_jList[[i]]$x_sq %*% theta$fixed_coef_hat  +
      t(theta_jList[[i]]$random_coef) %*% theta_jList[[i]]$z_sq %*% theta_jList[[i]]$random_coef -
      2 * theta_jList[[i]]$xy %*% theta$fixed_coef_hat  -
      2 * t(theta_jList[[i]]$zy) %*%  theta_jList[[i]]$random_coef +
      2 * t(theta$fixed_coef_hat) %*% theta_jList[[i]]$zx_mat %*%  theta_jList[[i]]$random_coef +
      theta$resid_var * sum(diag(theta_jList[[i]]$c_inv %*% theta_jList[[i]]$z_sq))

    theta$t1                     <- theta$t1 + theta_jList[[i]]$t1_j
    theta$t2                     <- theta$t2 + theta_jList[[i]]$t2_j
    theta$t3                     <- theta$t3 + theta_jList[[i]]$t3_j

  }
  theta$fixed_coef_hat <-
    theta$x_inv %*% (theta$xy_vector - theta$t1)
  theta$random_var_hat <- theta$t2 / J
  theta$resid_var_hat  <- theta$t3 / n
  class(theta) <- c("list", "sema")
  results <- list(model = theta,
                 unit = theta_jList)
 class(results) <- c("list", "sema")
 return(results)
}

predictY <- function(data.fixed,
                     theta_j,
                     theta,
                     data.random){
  if(is.null(theta_j)){
    theta_j <- list(random_coef = rep(0, length(data.random)))
  }
  return(
    as.matrix(data.fixed) %*% theta$fixed_coef_hat +
      as.matrix(data.random) %*% theta_j$random_coef
  )
}

lmerPredict	<- function(fixed_coef,
                        random_coef,
                        data.fixed,
                        data.random){
  if(dim(random_coef)[1]==0 )
  {
    random_coef <- matrix(0, ncol=length(data.random))
  }
  return(as.matrix(data.fixed) %*% fixed_coef +
          as.matrix(data.random) %*% as.matrix(t(random_coef)))
}

returnCovar <- function(covarMatrix){
  if(nrow(covarMatrix) == 1){
    return(covarMatrix)
  }
  else{
    return(c(diag(covarMatrix),
                 covarMatrix[lower.tri(covarMatrix)]))
  }
}
