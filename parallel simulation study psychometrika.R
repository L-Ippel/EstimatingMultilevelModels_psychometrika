################################################################################
#### SEMA for multilevel models
#### including random intercepts and slopes
#### and level 1 effects
#### psychometrika - round 3 2018
#### simulation study parallel
################################################################################
path <- ""
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
  path = args[1]
}

library(parallel)
library(SEMA)

#pre study definitions


### Conditions of simulations
iter = 1000
correlationsRandomVar <- rep(c(NA, 0, .15, .5), iter)
conditionsSim	<- cbind(correlationsRandomVar, iter = rep(1:iter, each = 4))
ncores		<- detectCores() / 2 #or the number i can use 

do.sim	<- function(pos, conditions = conditionsSim, path=path)
{
  source("simulation functions.R")
  set.seed(27840*pos)
  
  ### Select condition
  cor.random <- conditions[pos, 1] #selecteer correlatie 
  if(is.na(cor.random)){ 
    data.random <- 3 #locatie in data frame waar random effects staan 
    random.var  <- 50 #true waarde variantie random effect
  }
  else{
    data.random <- 3:7 #locatie in data frame waar random effects staan 
    random.var <- c(50, 0.2, 0.6, 1.8, 5)  #true waarde variantie random effect
  }
  
  trainingsize <- 2000
  lengthStream <- 50000
  n.individual <- lengthStream/50
  tempRes      <- tempResUpdate <- list() #objects nodig voor sema/sema update
  evaluate     <- 1000 #hoe vaak EM algorithme en en sema update gaan updaten 
  store        <- 100 # hoe vaak de resultaten opslaan 
  Window       <- 10000 #hoeveel data punten SWEM mee neemt
  
  ### objects to store the results
  predictionsSema <-
    predictionsUpdate <-
    predictionsEM <-
    predictionsSWEM <- matrix(nrow = lengthStream - trainingsize, ncol = 3)
  
  fixedSema <-
    fixedUpdate <-
    fixedEM <-
    fixedSWEM <- matrix(nrow = (lengthStream - trainingsize)/store, ncol = 15)
  
  randomSema <-
    randomUpdate <-
    randomEM <-
    randomSWEM <- matrix(nrow = (lengthStream - trainingsize)/store,
                         ncol = length(data.random) +
                           (length(data.random)*(length(data.random) - 1)/2))
  resVarSema <-
    resVarUpdate <-
    resVarEM <-
    resVarSWEM <- rep(NA, (lengthStream - trainingsize)/store)
  
  #generate data
  simData <- genData(n              = lengthStream,
                     j              = n.individual,
                     n.level2       = 7,
                     n.level1       = 8,
                     fixed.coef     = c(100, seq(0.1, 5.3, .4)),
                     random.var     = random.var,
                     cor.random.var = cor.random,
                     resid.var      = 5)
  
  #training set 
  fullEM <- emAlgorithm(data         = simData[1:trainingsize, ],
                        id           = 1,
                        y            = 2,
                        start.fixed  = c(100, seq(0.1, 5.3, .4)),
                        start.random = random.var,
                        data.random  = data.random,
                        start.cor    = cor.random,
                        start.res    = 5,
                        max.iter     = 800,
                        crit.value   = .0001)
  #initialiseer objecten voor SEMA/Update/SWEM
  id.parameters   <- id.parametersUpdate <- fullEM$unit
  tempRes$model   <- tempResUpdate$model <- fullEM$model
  id.list         <- SWEM.id.list        <- fullEM$id_list
  SWEM            <- fullEM
  
  for(n in (1 + trainingsize):nrow(simData)){
    id       <- simData$id[n]
    idTag    <- which(id.list == id) #id in data != idTag
    SWid     <- which(SWEM.id.list == id) #aangezien SWEM maar deel gebruikt, gebruikt SWEM eigen id tags
    
    #retrieve records
    idRecord       <- try.record(id = idTag, id.param = id.parameters)
    idRecordUpdate <- try.record(id = idTag, id.param = id.parametersUpdate)
    idRecordEM     <- try.record(id = idTag, id.param = fullEM$unit)
    idRecordSWEM   <- try.record(id = SWid,  id.param = SWEM$unit)
    
    if(is.null(idRecordSWEM)){ #make records in case it doesnt exist 
      SWEM.id.list <- c(SWEM.id.list, id)
      SWid         <- length(SWEM.id.list)
      
      if(is.null(idRecord)){
        id.list <- c(id.list, id)
        idTag   <- length(id.list)
      }
    }
    #predictions
    predictionsSema[n - trainingsize, ]   <- 
      c(id, simData$y[n], predictY(data.fixed  = simData[n, -c(1, 2)],
                                   theta_j     = idRecord,
                                   theta       = tempRes$model,
                                   data.random = simData[n, data.random]))
    predictionsUpdate[n - trainingsize, ] <- 
      c(id, simData$y[n], predictY(data.fixed  = simData[n, -c(1, 2)],
                                   theta_j     = idRecordUpdate,
                                   theta       = tempResUpdate$model,
                                   data.random = simData[n, data.random]))
    predictionsEM[n - trainingsize, ]     <- 
      c(id, simData$y[n], predictY(data.fixed  = simData[n, -c(1, 2)],
                                   theta_j     = idRecordEM,
                                   theta       = fullEM$model,
                                   data.random = simData[n, data.random]))
    predictionsSWEM[n - trainingsize, ]   <- 
      c(id, simData$y[n], predictY(data.fixed  = simData[n, -c(1, 2)],
                                   theta_j     = idRecordSWEM,
                                   theta       = SWEM$model,
                                   data.random = simData[n, data.random]))
    
    if((n %% evaluate) == 0){
      cat(c(n, conditions[pos, ]), '\n')
      
      fullUpdate <- semaUpdate(theta_jList = id.parametersUpdate,
                               theta       = tempResUpdate$model,
                               J           = tempResUpdate$model$j,
                               n           = tempResUpdate$model$n )
      tempResUpdate$model <- fullUpdate$model
      id.parametersUpdate <- fullUpdate$unit
      
      if(n != lengthStream){
        fullEM <- emAlgorithm(data         = simData[1:n, ],
                              id           = 1,
                              y            = 2,
                              start.fixed  = fullEM$model$fixed_coef_hat,
                              start.random = fullEM$model$random_var_hat,
                              data.random  = data.random,
                              start.cor    = cor.random,
                              start.res    = fullEM$model$resid_var_hat,
                              max.iter     = 20,
                              crit.value   = .0001)
      }
      if( n == lengthStream){ #extra iteraties aan t einde van de stroom 
        fullEM <- emAlgorithm(data         = simData[1:n,],
                              id           = 1,
                              y            = 2,
                              start.fixed  = fullEM$model$fixed_coef_hat,
                              start.random = fullEM$model$random_var_hat,
                              data.random  = data.random,
                              start.cor    = cor.random,
                              start.res    = fullEM$model$resid_var_hat,
                              max.iter     = 800,
                              crit.value   = .0001)
      }
      if(n <= Window){
        SWEM <- fullEM
      }
      else{
        SWEM <- emAlgorithm(data           = simData[(n - Window):n, ],
                            id             = 1,
                            y              = 2,
                            start.fixed    = SWEM$model$fixed_coef_hat,
                            start.random   = SWEM$model$random_var_hat,
                            data.random    = data.random,
                            start.cor      = cor.random,
                            start.res      = SWEM$model$resid_var_hat,
                            max.iter       = 20,
                            crit.value     = .0001)
      }
      SWEM.id.list <- SWEM$id_list
    }
    #SEMA
    tempRes <- sema_fit_one(data_fixed  = as.numeric(simData[n, -c(1, 2)]),
                            data_random = as.numeric(simData[n, data.random]),
                            data_y      = simData$y[n],
                            id          = id,
                            theta_j     = idRecord,
                            theta       = tempRes$model)
    
    #record terug in global SEMA list
    id.parameters[[idTag]] <- tempRes$unit
    #SEMA update
    tempResUpdate <- sema_fit_one(data_fixed  = as.numeric(simData[n, -c(1, 2)]),
                                  data_random = as.numeric(simData[n, data.random]),
                                  data_y      = simData$y[n],
                                  id          = id,
                                  theta_j     = idRecordUpdate,
                                  theta       = tempResUpdate$model)
    
    id.parametersUpdate[[idTag]] <- tempResUpdate$unit
    
    if((n %% store) == 0){ #store results 
      fixedSema[(n - trainingsize)/store, ]    <- tempRes$model$fixed_coef_hat
      fixedUpdate[(n - trainingsize)/store, ]  <- tempResUpdate$model$fixed_coef_hat
      fixedEM[(n - trainingsize)/store, ]      <- fullEM$model$fixed_coef_hat
      fixedSWEM[(n - trainingsize)/store, ]    <- SWEM$model$fixed_coef_hat
      
      randomSema[(n - trainingsize)/store, ]   <- returnCovar(tempRes$model$random_var_hat)
      randomUpdate[(n - trainingsize)/store, ] <- returnCovar(tempResUpdate$model$random_var_hat)
      randomEM[(n - trainingsize)/store, ]     <- returnCovar(fullEM$model$random_var_hat)
      randomSWEM[(n - trainingsize)/store, ]   <- returnCovar(SWEM$model$random_var_hat)
      
      resVarSema[(n - trainingsize)/store]     <- tempRes$model$resid_var_hat
      resVarUpdate[(n - trainingsize)/store]   <- tempResUpdate$model$resid_var_hat
      resVarEM[(n - trainingsize)/store]       <- fullEM$model$resid_var_hat
      resVarSWEM[(n - trainingsize)/store]     <- SWEM$model$resid_var_hat
      
    }
    
  }
  #save results
  write.table(predictionsSema, file = paste(path, "predSEMA_rv_", length(data.random),
                                            "_cor_", cor.random, "_it_", 
                                            conditions[pos, 2], ".txt", 
                                            sep = ""), sep = "   ")
  write.table(predictionsUpdate, file = paste(path, "predSU_rv_", length(data.random),
                                              "_cor_", cor.random, "_it_", 
                                              conditions[pos, 2], ".txt"
                                              , sep = ""), sep = "   ")
  write.table(predictionsEM, file = paste(path, "predEM_rv_", length(data.random),
                                          "_cor_", cor.random, "_it_", 
                                          conditions[pos, 2], ".txt",
                                          sep = ""), sep = "   ")
  write.table(predictionsSWEM, file = paste(path, "predSWEM_rv_", length(data.random),
                                            "_cor_", cor.random, "_it_", 
                                            conditions[pos, 2], ".txt", 
                                            sep = ""), sep = "   ")
  write.table(fixedSema, file = paste(path, "fixedSEMA_rv_", length(data.random),
                                      "_cor_", cor.random, "_it_", 
                                      conditions[pos, 2], ".txt", sep = ""), 
              sep = "   ")
  write.table(fixedUpdate, file = paste(path, "fixedSU_rv_", length(data.random),
                                        "_cor_", cor.random, "_it_", 
                                        conditions[pos, 2], ".txt"
                                        , sep = ""), sep = "   ")
  write.table(fixedEM, file = paste(path, "fixedEM_rv_", length(data.random),
                                    "_cor_", cor.random, "_it_", 
                                    conditions[pos, 2], ".txt"
                                    , sep = ""), sep = "   ")
  write.table(fixedSWEM, file = paste(path, "fixedSWEM_rv_", length(data.random),
                                      "_cor_", cor.random, "_it_", 
                                      conditions[pos, 2], ".txt"
                                      , sep = ""), sep = "   ")
  write.table(randomSema, file = paste(path, "randomSEMA_rv_", length(data.random),
                                       "_cor_", cor.random, "_it_", 
                                       conditions[pos, 2], ".txt"
                                       , sep = ""), sep = "   ")
  write.table(randomUpdate, file = paste(path, "randomSU_rv_", length(data.random),
                                         "_cor_", cor.random, "_it_", 
                                         conditionsSim[pos, 2], ".txt"
                                         , sep = ""), sep = "   ")
  write.table(randomEM, file = paste(path, "randomEM_rv_", length(data.random),
                                     "_cor_", cor.random, "_it_", 
                                     conditionsSim[pos, 2], ".txt"
                                     , sep = ""), sep = "   ")
  write.table(randomSWEM, file = paste(path, "randomSWEM_rv_", length(data.random),
                                       "_cor_", cor.random, "_it_", 
                                       conditions[pos, 2], ".txt"
                                       , sep = ""), sep = "   ")
  write.table(resVarSema, file = paste(path, "resVarSEMA_rv_", length(data.random),
                                       "_cor_", cor.random, "_it_", 
                                       conditions[pos, 2], ".txt"
                                       , sep = ""), sep = "   ")
  write.table(resVarUpdate, file = paste(path, "resVarSU_rv_", length(data.random),
                                         "_cor_", cor.random, "_it_", 
                                         conditions[pos, 2], ".txt"
                                         , sep = ""), sep = "   ")
  write.table(resVarEM, file = paste(path, "resVarEM_rv_", length(data.random),
                                     "_cor_", cor.random, "_it_", 
                                     conditions[pos, 2], ".txt"
                                     , sep = ""), sep = "   ")
  write.table(resVarSWEM, file = paste(path, "resVarSWEM_rv_", length(data.random),
                                       "_cor_", cor.random, "_it_", 
                                       conditions[pos, 2], ".txt"
                                       , sep = ""), sep = "   ")
  
  return(cat(conditions[pos,], '\n'))
}

if (ncores == 1){
  ### Sequential version of running simulation
  out <- lapply(1:nrow(conditionsSim), do.sim, conditions = conditionsSim)
} else{
  cl <- makeCluster(ncores, outfile="") # Create cluster
  clusterCall(cl, function() library(SEMA, mvtnorm))
  
  clusterExport(cl, varlist = c("conditionsSim"))
  out <- clusterApplyLB(cl, 1:nrow(conditionsSim), do.sim,
                        conditions = conditionsSim, path = path) # Run simulation
  stopCluster(cl) # Shut down the nodes
}
