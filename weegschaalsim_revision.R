#########################################
##scale data - replay the data stream - september 2018
#########################################
#load the data
myData                <- read.table("simWeightdata.txt", header = T, sep = ";")
myData$ageCentered    <- 1970 - myData$dateOfBirth
myData$lengthCentered <- myData$height - 174
myData$id             <- myData$nomem_encr
myData$y              <- myData$weight
myData                <- myData[, c(24, 25, 3:23)]
#Bevolking; kerncijfers, 03 oktober 2016 statline

#load SEMA function
source("simulation functions.R")


fixedVar <- c(3, 6, 22, 23, 7:19) 
randomVar <- c(3, 11:16)

#create objects to store results

# starting values Beta -> statline: Lengte en gewicht van personen, ondergewicht en overgewicht; vanaf 1981
# selection Opgegeven lengte, 20 jaar of ouder Opgegeven gewicht, 20 jaar of ouder

library(SEMA)
library(lme4)

m1 <- lmer(y ~ 1 + female + ageCentered+ lengthCentered + daily + 
             weekly + norm + target + sunday + monday + tuesday + wednesday + thursday + 
             saturday + night + afternoon + evening + 
             (1 +  sunday + monday + tuesday + wednesday + thursday + 
                saturday | id), data = myData[1:6894,], REML = FALSE)
#pre study definitions

trainingsize <- 6894
lengthStream <- nrow(myData)
tempRes      <- tempResUpdate <- list()
store        <- 100
Window       <- 12000

fullEM  <- emAlgorithm(data         = myData[1:trainingsize, c(1,2,fixedVar)],
                       id = 1, 
                       y = 2,
                       start.fixed   = as.numeric(fixef(m1)),
                       start.random  = as.data.frame(VarCorr(m1))$vcov[1:length(randomVar)], 
                       data.random   = randomVar,
                       start.cor     = 0.2, 
                       start.res     = 1,
                       max.iter      = 2000, 
                       crit.value    = .0001)

id.parameters   <- id.parametersUpdate <- fullEM$unit
tempRes$model   <- tempResUpdate$model <- fullEM$model
id.list         <- SWEM.id.list        <- fullEM$id_list
SWEM            <- fullEM

#objects to store the results
predictionsSema <- 
  predictionsUpdate <- 
  predictionsEM <- 
  predictionsSWEM <- matrix(nrow = lengthStream - trainingsize, ncol = 3)

fixedSema <- 
  fixedUpdate <- 
  fixedEM <- 
  fixedSWEM <- matrix(nrow = (lengthStream - trainingsize)/store, 
                      ncol = length(fixedVar))

randomSema <- 
  randomUpdate <- 
  randomEM <- 
  randomSWEM <- matrix(nrow = (lengthStream - trainingsize)/store, 
                       ncol = length(randomVar) + 
                         (length(randomVar)*(length(randomVar) - 1)/2))
resVarSema <- 
  resVarUpdate <- 
  resVarEM <- 
  resVarSWEM <- rep(NA, (lengthStream - trainingsize)/store)

for(n in (1 + 163000):nrow(myData)){
  id       <- myData$id[n]
  idTag    <- which(id.list == id) 
  SWid     <- which(SWEM.id.list == id)
  
  idRecord       <- try.record(id = idTag, id.param = id.parameters)
  idRecordUpdate <- try.record(id = idTag, id.param = id.parametersUpdate)
  idRecordEM     <- try.record(id = idTag, id.param = fullEM$unit)
  idRecordSWEM   <- try.record(id = SWid,  id.param = SWEM$unit)
  
  if(is.null(idRecordSWEM)){ 
    SWEM.id.list <- c(SWEM.id.list, id) 
    SWid         <- length(SWEM.id.list)
    
    if(is.null(idRecord)){
      id.list <- c(id.list, id)
      idTag   <- length(id.list)
    }
  }
  
  predictionsSema[n - trainingsize, ] <- c(id, 
                                           myData$y[n],
                                           predictY(theta = tempRes$model, 
                                                   theta_j     = idRecord,
                                                   data.fixed  = (myData[n, fixedVar]),
                                                   data.random = (myData[n, randomVar])))
  predictionsUpdate[n - trainingsize, ] <- c(id, 
                                             myData$y[n],
                                             predictY(theta = tempResUpdate$model, 
                                                     theta_j     = idRecordUpdate,
                                                     data.fixed   = (myData[n, fixedVar]),
                                                     data.random = (myData[n, randomVar])))
                                             
  predictionsEM[n - trainingsize, ]     <- c(id, 
                                             myData$y[n],
                                             predictY(theta = fullEM$model, 
                                                     theta_j     = idRecordEM,
                                                     data.fixed   = (myData[n, fixedVar]),
                                                     data.random = (myData[n, randomVar]))) 
                                            
  predictionsSWEM[n - trainingsize, ]   <- c(id, 
                                             myData$y[n],
                                             predictY(theta = SWEM$model, 
                                                     theta_j     = idRecordSWEM,
                                                     data.fixed   = (myData[n, fixedVar]),
                                                     data.random = (myData[n, randomVar]))) 
  if( n == 163000){ #retrain the model after new units entered 
    fullEM <- emAlgorithm(data         = myData[1:n, c(1,2,fixedVar)],
                          id = 1, 
                          y = 2,
                          start.fixed  = fullEM$model$fixed_coef_hat,
                          start.random = fullEM$model$random_var_hat, 
                          data.random  = randomVar,
                          start.cor    = .2, 
                          start.res    = fullEM$model$resid_var_hat,
                          max.iter     = 2000, 
                          crit.value   = .0001)
    id.parameters   <- id.parametersUpdate <- fullEM$unit
    tempRes$model   <- tempResUpdate$model <- fullEM$model
    id.list         <- SWEM.id.list        <- fullEM$id_list
    SWEM            <- fullEM
    
  } 
  
  if(n>(trainingsize+1)& (myData$weighdate[n] !=  myData$weighdate[n-1])){
    for(k in 1:2){#extra iteration to learn a bit more
      fullUpdate <- semaUpdate(theta_jList  = id.parametersUpdate, 
                               theta       = tempResUpdate$model,
                               J           = tempResUpdate$model$j, 
                               n           = tempResUpdate$model$n )
      tempResUpdate$model <- fullUpdate$model
      id.parametersUpdate <- fullUpdate$unit
    } 
    
    if(n <= Window){
      SWEM <- fullEM
    }
    else{
      SWEM <- emAlgorithm(data           = myData[(n - Window):n, c(1,2,fixedVar)],
                          id = 1, 
                          y = 2,
                          start.fixed    = SWEM$model$fixed_coef_hat,
                          start.random   = SWEM$model$random_var_hat, 
                          data.random    = randomVar,
                          start.cor      = .2, 
                          start.res      = SWEM$model$resid_var_hat,
                          max.iter       = 100, 
                          crit.value     = .0001)
    }
    SWEM.id.list <- SWEM$id_list
    
  }
  
     
  if(myData$sunday[n] == 1 & myData$sunday[n+1] == 0 ){  
      fullEM <- emAlgorithm(data         = myData[1:n, c(1,2,fixedVar)],
                            id = 1, 
                            y = 2,
                            start.fixed  = fullEM$model$fixed_coef_hat,
                            start.random = fullEM$model$random_var_hat, 
                            data.random  = randomVar,
                            start.cor    = .2, 
                            start.res    = fullEM$model$resid_var_hat,
                            max.iter     = 1000, 
                            crit.value   = .0001)
    
      
    }
  
  tempRes <- sema_fit_one(data_fixed  = as.numeric(myData[n, fixedVar]),
                          data_random = as.numeric(myData[n, randomVar]),
                          data_y      = myData$y[n],
                          id          = id,
                          theta_j     = idRecord, 
                          theta       = tempRes$model)
  
  
  id.parameters[[idTag]] <- tempRes$unit
  
  tempResUpdate <- sema_fit_one(data_fixed  = as.numeric(myData[n, fixedVar]),
                                data_random = as.numeric(myData[n, randomVar]),
                                data_y      = myData$y[n],
                                id          = id,
                                theta_j     = idRecordUpdate, 
                                theta       = tempResUpdate$model)
  
  id.parametersUpdate[[idTag]] <- tempResUpdate$unit
  
  if((n %% store) == 0){
    print(n)
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
write.table(predictionsSema,   file = "ScalepredSEMAr2.txt",   sep = "   ")
write.table(predictionsUpdate, file = "ScalepredSUr2.txt",     sep = "   ")
write.table(predictionsEM,     file = "ScalepredEMr2.txt",     sep = "   ")
write.table(predictionsSWEM,   file = "ScalepredSWEMr2.txt",   sep = "   ")
write.table(fixedSema,         file = "ScalefixedSEMAr2.txt",  sep = "   ")
write.table(fixedUpdate,       file = "ScalefixedSUr2.txt",    sep = "   ")
write.table(fixedEM,           file = "ScalefixedEMr2.txt",    sep = "   ")
write.table(fixedSWEM,         file = "ScalefixedSWEMr2.txt",  sep = "   ")
write.table(randomSema,        file = "ScalerandomSEMAr2.txt", sep = "   ")
write.table(randomUpdate,      file = "ScalerandomSUr2.txt",   sep = "   ")
write.table(randomEM,          file = "ScalerandomEMr2.txt",   sep = "   ")
write.table(randomSWEM,        file = "ScalerandomSWEMr2.txt", sep = "   ")
write.table(resVarSema,        file = "ScaleresVarSEMAr2.txt", sep = "   ")
write.table(resVarUpdate,      file = "ScaleresVarSUr2.txt",   sep = "   ")
write.table(resVarEM,          file = "ScaleresVarEMr2.txt",   sep = "   ")
write.table(resVarSWEM,        file = "ScaleresVarSWEMr2.txt", sep = "   ")

