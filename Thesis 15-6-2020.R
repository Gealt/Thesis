rm(list = ls())
library(doParallel)
library(foreach)

# #one event ID per specific dyad
# #statlist = rows -> eventtimes, columns -> dyads, slice -> variables(intercept = all 1's, variables)

#array[3][200][events][pairs][sufficientstatistics]
# baseline, 200 datasets, eventtime, pairs, sufficientstatistics 

set.seed(14-4-2020)
updateEndogenous <- function(tuple,pairs){
  Sender <- tuple[nrow(tuple),"Sender"]
  Receiver <- tuple[nrow(tuple), "Receiver"]
  pair <- paste0(Sender,",",Receiver)
  
  #all pairs with the same sender, excluding the current pair
  temp <- pairs[pairs != pair][grepl(paste0(Sender, ","), pairs[pairs != pair])]
  
  endogenous[pair, "Persistence"] <- (sum(tuple[,"Sender"] == Sender & tuple[,"Receiver"] == Receiver)
                                      /sum(tuple[,"Sender"] == Sender))
  endogenous[temp, "Persistence"] <- (endogenous[temp, "Persistence"] * (sum(tuple[,"Sender"] == Sender)-1)
                                      /sum(tuple[,"Sender"] == Sender))
  
  endogenous[pair, "Preferential Attachment"] <- (sum(tuple[,"Sender"] == Receiver | tuple[,"Receiver"] == Receiver)
                                                  /(nrow(tuple)*2))
  
  endogenous[pair, "Recency"] <- 1
  endogenous[temp, "Recency"] <- 1/(order(1/endogenous[temp, "Recency"])+1)
  
  BA <- paste0(Receiver,",",Sender)
  endogenous[, "AB-BA"] <- ifelse(names(endogenous[,"AB-BA"]) == BA,1,0)
  
  BY <- pairs[pairs != paste0(Receiver, ",", Sender)][grepl(paste0(Receiver, ","), 
                                                            pairs[pairs != paste0(Receiver, ",", Sender)])]
  endogenous[, "AB-BY"] <- ifelse(names(endogenous[,"AB-BY"]) %in% BY,1,0)
  
  XA <- pairs[pairs != paste0(Receiver, ",", Sender)][grepl(paste0(",", Sender), 
                                                            pairs[pairs != paste0(Receiver, ",", Sender)])]
  endogenous[, "AB-XA"] <- ifelse(names(endogenous[,"AB-XA"]) %in% XA,1,0)
  
  XB <- pairs[pairs != pair][grepl(paste0(",", Receiver), 
                                   pairs[pairs != pair])]
  endogenous[, "AB-XB"] <- ifelse(names(endogenous[,"AB-XB"]) %in% XB,1,0)
  
  AB <- pair
  endogenous[, "AB-AB"] <- ifelse(names(endogenous[,"AB-AB"]) == AB,1,0)
  
  AY <- pairs[pairs != pair][grepl(paste0(Sender, ","), 
                                   pairs[pairs != pair])]
  endogenous[, "AB-AY"] <- ifelse(names(endogenous[,"AB-AY"]) %in% AY,1,0)
  
  XY <- endogenous[, c("AB-BA", "AB-BY", "AB-XA", "AB-XB",
                       "AB-AB", "AB-AY")] == 0
  endogenous[, "AB-XY"] <- apply(XY,1, all)
  
  
  #add the triadic effects for the non-current pair aswell
  for (Sender in names){
    for (Receiver in names){
      if (Sender == Receiver){
        next
      }
      pair <- paste0(Sender,",",Receiver)
      sumO2P <- 0 #sum outgoing 2 paths
      sumI2P <- 0 #sum incoming 2 paths
      sumOSP <- 0 #sum outgoing shared partners
      sumISP <- 0 #inbound shared partners
      for (possibleReceiver in names){
        sumO2P <- sumO2P+min(sum(tuple[,"Sender"] == Sender & tuple[,"Receiver"] == possibleReceiver),
                             sum(tuple[,"Sender"] == possibleReceiver & tuple[,"Receiver"] == Receiver)) 
        sumI2P <- sumI2P+min(sum(tuple[,"Sender"] == Receiver & tuple[,"Receiver"] == possibleReceiver),
                             sum(tuple[,"Sender"] == possibleReceiver & tuple[,"Receiver"] == Sender))
        sumOSP <- sumOSP+min(sum(tuple[,"Sender"] == Sender & tuple[,"Receiver"] == possibleReceiver),
                             sum(tuple[,"Sender"] == Receiver & tuple[,"Receiver"] == possibleReceiver))
        sumISP <- sumISP+min(sum(tuple[,"Sender"] == possibleReceiver & tuple[,"Receiver"] == Sender),
                             sum(tuple[,"Sender"] == possibleReceiver & tuple[,"Receiver"] == Receiver))
      } 
      endogenous[pair, "Outgoing2Paths"] <- sumO2P
      endogenous[pair, "Incoming2Paths"] <- sumI2P
      endogenous[pair, "OutgoingSharedPartner"] <- sumOSP
      endogenous[pair, "InboundSharedPartner"]  <- sumISP 
    }
  }
  

  endogenous[, "Outgoing2Paths"] <- endogenous[, "Outgoing2Paths"]/max(endogenous[, "Outgoing2Paths"])
  endogenous[, "Incoming2Paths"] <- endogenous[, "Incoming2Paths"]/max(endogenous[, "Incoming2Paths"])
  endogenous[, "OutgoingSharedPartner"] <- endogenous[, "OutgoingSharedPartner"]/max(endogenous[, "OutgoingSharedPartner"])
  endogenous[, "InboundSharedPartner"]  <- endogenous[, "InboundSharedPartner"]/max(endogenous[, "InboundSharedPartner"])
  endogenous[is.nan(endogenous)] <- 0
    
  return (endogenous)
}

getpairs <- function(names){
  #initialize an empty vector to hold all possible combinations of pairs
  pairs <-vector()
  #fill the pairs vector 
  #excluding pairs in which a participant is both sender and receiver
  for (sender in names){
    for (receiver in names){
      if (sender != receiver)
        pairs <- c(pairs, paste0(sender,",",receiver))
    }
  }
  return(pairs)
}

initialiseExogenous <- function(names){
  #initialise a matrix for the exogenous variables with random assigned roles
  #e.g. participants ones roles are randomly assigned to be a working person who follows education
  exogenous <- matrix(sample(c(0,1), length(names)*3, TRUE),
                      nrow = length(names), ncol = 3)
  rownames(exogenous) <- names
  colnames(exogenous) <- c("work", "education", "sport")
  return(exogenous)
}

initialiseEndogenous <- function(names,pairs){
  #initialise a matrix for the exogenous variables, initialised with the intercept and 0 starting values
  endogenous <- matrix(0, nrow = length(pairs), ncol = 15)
  endogenous[,1] <- 1
  rownames(endogenous) <- pairs
  colnames(endogenous) <- c("Intercept","Persistence", "Preferential Attachment", "Recency", 
                            "Outgoing2Paths", "Incoming2Paths", "OutgoingSharedPartner",
                            "InboundSharedPartner", "AB-BA", "AB-BY", "AB-XA",
                            "AB-XB", "AB-XY", "AB-AB", "AB-AY")
  return(endogenous)
}


cores=detectCores() #check how many cores the computer has
cl <- makeCluster(cores[1]-1) #use all but one core to prevent overloading of cpu
registerDoParallel(cl) #register the cluster for use

#split 3 baselines in 6 components, start with the longest event sequences <- to_do

#loop over each baseline in parallel, combining them by their columns
names <- c(1, 2, 3, 4, 5, 6, 7)
pairs <- getpairs(names)
exogenous <- initialiseExogenous(names)
test <- foreach(sequencelength = c(750,750,750,500,500,500,250,250,250), .combine = cbind, .packages = c("abind")) %dopar%{
  baseline <- -5
  #storage vector for all datasets in the current baseline
  allDataSets <- vector(67, mode = "list")
  
  for (i in 1:length(allDataSets)){
    #Possible senders / receivers
    endogenous <- initialiseEndogenous(names,pairs)
    
    #initialise sufficientstatistics for the current event
    eventSufficientStatistics <- array(NA,c(length(pairs),21))
    rownames(eventSufficientStatistics) <- pairs
    colnames(eventSufficientStatistics) <- c(paste0(colnames(exogenous),"Sender"),paste0(colnames(exogenous),"Receiver"),colnames(endogenous))
    
    #initialise the eventlist and statlist
    eventlist <- c(0,0)
    statlist <- array(0,c(1,dim(eventSufficientStatistics)))
    dimnames(statlist) <- list(0,rownames(eventSufficientStatistics),colnames(eventSufficientStatistics))
    
    for (pair in pairs){
      sender <- strsplit(pair, ",")[[1]][1]
      receiver <- strsplit(pair, ",")[[1]][2]
      statlist[1,pair,c(paste0(colnames(exogenous),"Sender"),paste0(colnames(exogenous),"Receiver"))] <- c(exogenous[sender,],exogenous[receiver,])
    }    
    
    statlist[1,,"Intercept"] <- 1
    #Thetas for the sufficientstatistics (for when generating the true model
    #make sure to put several thetas at 0)
    theta <- c(-1.4,0.7,0,  #exogenous variables of the sender
               0,0,0,   #exogenous variables of the receiver
               0,0,1.48,1.34, #"Intercept", Persistence", "Preferential Attachment", "Recency"
               0.0,0.7,0,    #"Outgoing2Paths", "Incoming2Paths", "OutgoingSharedPartner"
               0, 0, 1.62, 0, #"InboundSharedPartner", "AB-BA", "AB-BY", "AB-XA",
               0, 1.75, 0, 0) #"AB-XB", "AB-XY", "AB-AB", "AB-AY"
    
    #AB-XY, O2P, Persistence, Recency, AB-BY,work,school
    
    #initialise a tuple of event time, sender and receiver --> remove the initial row 
    tuple <- data.frame(Time = 0, Sender = 1, Receiver = 2)
    
    #initialise empty vector for the lograte
    lograte<- vector()
    
    while (nrow(tuple) < sequencelength){
      #update endogenous at the start, because the tuple is initialised with an event
      endogenous <- updateEndogenous(tuple,pairs)
      
      
      for (pair in pairs){
        #get the sender/receiver out of the current pair
        sender <- strsplit(pair, ",")[[1]][1]
        receiver <- strsplit(pair, ",")[[1]][2]
        #make sufficientstatistics for the current pair
        pairSufficientStatistics <- c(exogenous[sender,],exogenous[receiver,],
                                      endogenous[pair,])
        eventSufficientStatistics[pair,] <- pairSufficientStatistics
        #fill log rates
        lograte[pair] <- baseline+sum(theta*eventSufficientStatistics[pair,])
      }
      #make the lograte log-linear
      rate <- exp(lograte)
      #sum the log linear logrates
      sumrate <- sum(rate)
      #sample an intertime time
      intertime <- rexp(n = 1, rate = sumrate)
      #Assign  time for the latest event tuple
      Time <- tuple[nrow(tuple),"Time"]+intertime
      #sample a pair from the possible pairs
      pair <- sample(pairs, size = 1, prob = rate/sumrate)
      #assign the current sender to event tuple
      Sender <- strsplit(pair, ",")[[1]][1]
      #assign the current receiver to event tuple
      Receiver <- strsplit(pair, ",")[[1]][2]
      #add latest event tuple to the total list of event tuples
      tuple <- rbind(tuple, data.frame(Time,Sender,Receiver))
      
      #add the current event to the eventlist, which is compatible with rem
      eventlist <- rbind(eventlist,c(match(pair,pairs),Time))
      
      #add the current eventsufficientstatistics to the
      statlist <- abind(statlist,eventSufficientStatistics, along = 1)
      dimnames(statlist)[[1]][nrow(statlist[,,1])] <- Time
    }
    
    #add the current dataset to the list of datasets for the current baseline
    allDataSets[[i]] <- list(eventlist,statlist)
  }
  allDataSets #cbinds the datasets in this baseline to the other baselines
  #to do: change variablename test
}

stopCluster(cl)
tmp <- c(0,0,0,0,0,0,0,0,0)
for (i in 1:ncol(test)){for (j in 1:nrow(test)){tmp[i] <- tmp[i] + nrow(test[[j,i]][[1]])}}


library(relevent)

#function to get effects for models in a bottom up fashion
GetSuperModels <- function(all_effects, old_models = vector(mode = "list")){
  if (!length(old_models)) return (all_effects) #if there are no models yet return effects
  #get all combinations of effects with the old models
  supermodel <- expand.grid(names(old_models), all_effects, stringsAsFactors = FALSE)
  #remove combinations with repeated effects
  supermodel <- supermodel[!(mapply(grepl,supermodel[,2],supermodel[,1], MoreArgs = list(fixed = TRUE))),]
  #take the combinations of effects into 1 string
  supermodel <- paste0(supermodel[,1],",",supermodel[,2])
  #remove double combinations in a different ordern
  tmp <- t(sapply(strsplit(supermodel,","),sort))
  for (row in 1:nrow(tmp)){
    string <- tmp[row,1]
    for (index in 2:ncol(tmp)){
      string <- paste0(string,",",tmp[row,index])
    }
    supermodel[row] <- string
  }
  return (unique(supermodel))
}

allModels <- function(array, all_effects, estimator = "MLE"){
  all_models <- vector(mode = "list")
  old_models <- vector(mode = "list")
  minBIC <- Inf
  while(tryCatch({
    GetSuperModels(all_effects, old_models)
    1
  }, 
  error = function(cond) return(0))){
    #reset the new models vector
    new_models <- vector(mode = "list")
    #run the rem function for each effect for the new models
    for (i in GetSuperModels(all_effects, old_models)){
      statlist <- array[[2]][,,c("Intercept",unlist(strsplit(i,",")))]
      new_models[[i]] <- rem(eventlist = array[[1]],
                             statslist = statlist,
                             estimator = estimator,
                             timing = "interval")
      if (new_models[[i]]$BIC < minBIC){
        minBIC <- new_models[[i]]$BIC
      }
    }
    old_models <- new_models
    all_models <- c(all_models,old_models)
  }
  for (i in names(all_models)){
    all_models[[i]][["marginal"]] <- exp(-(all_models[[i]]$BIC-minBIC)/2)
  }
  return (all_models)
}

validModels <- function(all_models, threshold = 20){
  maxposterior <- 0
  sumposterior <- 0
  prior <- 1/length(all_models)
  for (i in names(all_models)){
    posterior <- all_models[[i]]$marginal*prior
    sumposterior <- sumposterior + posterior
  }
  for (i in names(all_models)){
    posterior <- all_models[[i]]$marginal*prior/sumposterior
    if (posterior > maxposterior){
      maxposterior <- posterior
    }
  }
  for (i in names(all_models)){
    posterior <- all_models[[i]]$marginal*prior/sumposterior
    if (maxposterior/posterior > threshold){
      all_models[[i]] <- NULL
    }
  }
  return (all_models)
}

AverageModel <- function(models){
  prior <- 1/length(models)
  posterior <- vector(mode = "integer")
  averagemodel <- matrix(0, nrow = 2, ncol = length(all_effects)+2)
  colnames(averagemodel) <- c(all_effects,"Intercept", "Modelcount")
  for (i in names(models)){
    posterior[i] <- models[[i]]$marginal*prior
    coef <- models[[i]]$coef
    averagemodel[1,names(coef)] <- averagemodel[1,names(coef)]+coef*posterior[i]
  }
  averagemodel[1,] <- averagemodel[1,]/sum(posterior)
  averagemodel[1,"Modelcount"] <- length(models)
#  print(apply(bootstrapBill(models, posterior)$t, 2, mean))
  averagemodel[2,] <- apply(bootstrapBill(models, posterior)$t, 2, mean)
#  averagemodel[3,] <- apply(bootstrapWilliam(models,posterior)$t, 2, mean)
  return (averagemodel)
}

getSD <- function(data, index){
  averagemodel <- matrix(0, nrow = length(index), ncol = length(all_effects)+2)
  colnames(averagemodel) <- c(all_effects,"Intercept", "Modelcount")
  counter <- 0
  for (i in index){
    counter <- counter + 1
    averagemodel[counter,names(data[[i]]$coef)] <- data[[i]]$coef
  }
  return (apply(averagemodel, 2, sd))
}

getMean <- function(data,index){
  averagemodel <- matrix(0, nrow = length(index), ncol = length(all_effects)+2)
  colnames(averagemodel) <- c(all_effects,"Intercept", "Modelcount")
  counter <- 0
  for (i in index){
    counter <- counter + 1
    averagemodel[counter,names(data[[i]]$coef)] <- data[[i]]$coef
  }
  return (apply(averagemodel, 2, mean))
}

bootstrapBill <- function(models, posterior){
  boot(models, getSD, weights = posterior, R = 20000)
}

bootstrapWilliam <- function(models, posterior){
  boot(models, getMean, weights = posterior, R = 50000)
}

all_effects <- c("Incoming2Paths","AB-XY","AB-BY","AB-AY","Preferential Attachment", "Recency",
                 "workSender","educationSender")

for (estimator in c("BPM","MLE")){
  cl <- makeCluster(cores[1]-1) #use all but one core, to prevent CPU overloading
  registerDoParallel(cl) #register the cores for use
  models <- foreach(column = 1:ncol(test), .combine = cbind, .packages = c("relevent","boot")) %dopar%{
    for (row in 1:nrow(test)){
      set.seed(14-4-2020)
      allmodels <- allModels(test[[row,column]], all_effects, estimator = estimator)
      subset1 <- validModels(allmodels)
      all <- AverageModel(allmodels)
      occam1 <- AverageModel(subset1)
      if (column <= 3){
        if (estimator == "BPM") next
        saveRDS(all, file = paste0("C:/Users/Enrico/Desktop/thesis/",estimator,"/750/all/",row,"_",column))
        saveRDS(occam1, file = paste0("C:/Users/Enrico/Desktop/thesis/",estimator,"/750/occam/",row,"_",column))
      } else if (column > 3 & column <= 6){
        saveRDS(all, file = paste0("C:/Users/Enrico/Desktop/thesis/",estimator,"/500/all/",row,"_",column))
        saveRDS(occam1, file = paste0("C:/Users/Enrico/Desktop/thesis/",estimator,"/500/occam/",row,"_",column))
      } else {
        saveRDS(all, file = paste0("C:/Users/Enrico/Desktop/thesis/",estimator,"/250/all/",row,"_",column))
        saveRDS(occam1, file = paste0("C:/Users/Enrico/Desktop/thesis/",estimator,"/250/occam/",row,"_",column))
      }
    }
  }
  stopCluster(cl) #remove connection from cores
}

printstatistics <- function(basepath, statistics = "sd", estimator = "BSIR"){
  averages <- array(dim = c(3,3,2,201,2,length(all_effects)+2))
  dimnames(averages)[[1]] <- c("MLE","BPM","BSIR")
  dimnames(averages)[[2]] <- c(250,500,750)
  dimnames(averages)[[3]] <- c("all","occam")
  dimnames(averages)[[4]] <- c(1:201)
  dimnames(averages)[[5]] <- c("BMA","bootstrap")
  dimnames(averages)[[6]] <- c(all_effects,"Intercept","Modelcount")

  for (est in list.files(path=basepath)){
    for (sequence in list.files(paste0(basepath,"/",est))){
      for (BMAtype in list.files(paste0(basepath,"/",est,"/",sequence))){
        count <- 0
        for (file in list.files(paste0(basepath,"/",est,"/",sequence,"/",BMAtype))){
          count <- count + 1
          averages[est,sequence,BMAtype,count,,] <- readRDS(paste0(basepath,"/",est,"/",sequence,"/",BMAtype,"/",file))
        }
        if (est == estimator){
          print (c(estimator,sequence,BMAtype))
          if (statistics == "mean"){
            print(colMeans(averages[estimator,sequence,BMAtype,,1,]))
          }
          else if (statistics == "sd"){
            print(apply(averages[estimator,sequence,BMAtype,,1,], 2, sd))
          }
        }
      }
      if (statistics == "ttests"){
        print (c(estimator,sequence))
        print (ttests() < 0.05)
      }
    }
  }
  return (averages)
}

ttests <- function(){
  tmp <- vector(length(all_effects)+1, mode = "integer")
  names(tmp) <- colnames(averages[estimator,sequence,"occam",,])
  for (i in colnames(averages[estimator,sequence,"occam",,])){
    if (i == "Intercept") next
    tmp[i] <- t.test(averages[estimator,sequence,"all",,i], averages[estimator,sequence,"occam",,i])$p.value
  }
  return (tmp)
}

for (est in list.files(path=basepath)){
  for (sequence in list.files(paste0(basepath,"/",est))){
    for (BMAtype in list.files(paste0(basepath,"/",est,"/",sequence))){
      for (file in list.files(paste0(basepath,"/",est,"/",sequence,"/",BMAtype),full.names = TRUE)){
        print(file)
      }
    }
  }
}


test <- readRDS(file= "C:/Users/Enrico/Desktop/thesis/Ricardo")
all_effects <- c("Incoming2Paths","AB-XY","AB-BY","AB-AY","Preferential Attachment", "Recency",
                 "workSender","educationSender")

microbenchmark(rem(eventlist = test[[1]],
                                       statslist = test[[2]][,,c(all_effects,"Intercept")],
                                       estimator = "BSIR",
                                       timing = "interval"), times = 5)
