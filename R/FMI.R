#' On the conventional definition of path-specific effects - fully mediated interaction with multiple ordered mediators.
#'
#' FMI calculates the causal path-specific effects by decomposing fully mediated interaction.
#'
#' @name FMI
#' @author An-Shun Tai \email{anshuntai@nctu.edu.tw}
#' @param data A data frame, where the column is the variable and row is the sample.
#' @param exposure Variable name of the treatment or exposure.
#' @param mediators A vector of ordered mediator names.
#' @param outcome Variable name of the outcome.
#' @param baseline.conf A vector of confounder names.
#' @param med.model A vector of model names (GLM) to specify the mediator models.
#' @param outcome.model The outcome model
#' @return Estimated PSEs, FMI, and total effect (TE)
#' @export
#' @examples
#' FMI(data=sim.data,outcome="feeltherm",treat = "spifev1per",mediators=c("m1_crqdysmnP","m2_hadsanxN"),
#' baseline.conf = c("countryB","female","ageT0","bmi","smopackyrs"),
#' med.model=c("gaussian","gaussian"),
#' outcome.model = "gaussian",nboost=100)


FMI <- function(data,outcome = "outcome.name",treat = "treat.name",
                     mediators, baseline.conf = c("BaselineConf.name1","BaselineConf.name2",...),
                     med.model,
                     outcome.model,
                     a0 = 0,a1 = 1, MonteCarlo.time = 10^4,seed=1234,nboost=1000){
  FMI_main <- function(data,outcome,treat,
                       mediators, baseline.conf,
                       med.model,
                       outcome.model,
                       a0,a1, MonteCarlo.time){

    bconf.num <- length(baseline.conf)
    med.num <- length(mediators)

    # Step 1: model fitting
    for(t in 1:med.num ){
      ifelse(t==1,ind.med.t <- 0,ind.med.t <- 1:(t-1))

      # mediators
      form0 <- paste( c(mediators[ind.med.t],treat),collapse="*")
      form1 <- paste( c(form0,baseline.conf),collapse=" + ")
      form.med <- paste( c(mediators[t]," ~ ", form1), collapse="")
      formula <- as.formula(form.med)
      assign(paste("model_M",t,sep=""),glm(formula,family=med.model[t],data=data))
    }

    # outcome
    form0 <- paste( c( mediators, treat), collapse="*")
    form1 <- paste( c( form0, baseline.conf), collapse=" + ")
    form.outcome <- paste( c(outcome," ~ ", form1), collapse="")
    formula <- as.formula(form.outcome)
    modelY <- glm(formula,family=outcome.model ,data=data)


    # Second step
    phi0 <- phi1 <- rep(NA,1)
    phi <- rep(NA,med.num+1)

    data_MonteCarlo <- data[sample(1:nrow(data),MonteCarlo.time,replace =T),]
    data_MonteCarlo_s1 <- data_MonteCarlo_s0 <- data_MonteCarlo
    data_MonteCarlo_s1[,treat] <- a1
    for(t in 1:med.num ){
      # mediators
      data_MonteCarlo_s1[,mediators[t]] <- predict(get(paste("model_M",t,sep="")),data_MonteCarlo_s1,type = "response")
    }
    y_hat <- predict(modelY,data_MonteCarlo_s1,type = "response")
    phi1 <- mean(y_hat)

    data_MonteCarlo_s0[,treat] <- a0
    for(t in 1:med.num ){
      # mediators
      data_MonteCarlo_s0[,mediators[t]] <- predict(get(paste("model_M",t,sep="")),data_MonteCarlo_s0,type = "response")
    }
    y_hat <- predict(modelY,data_MonteCarlo_s0,type = "response")
    phi0 <- mean(y_hat)

    #
    data_MonteCarlo_temp <- data_MonteCarlo_s0
    data_MonteCarlo_temp[,treat] <- a1
    y_hat <- predict(modelY,data_MonteCarlo_temp,type = "response")
    phi[med.num+1] <- mean(y_hat)

    #
    for(p in 1: (med.num)){
      data_MonteCarlo_temp <- data_MonteCarlo_s0
      data_MonteCarlo_temp[,treat] <- a1

      data_MonteCarlo_temp[,mediators[p]] <- predict(get(paste("model_M",p,sep="")),data_MonteCarlo_temp,type = "response")
      if(p<med.num){
        for(ps in (p+1):(med.num)){
          data_MonteCarlo_temp[,treat] <- a0
          data_MonteCarlo_temp[,mediators[ps]] <- predict(get(paste("model_M",ps,sep="")),data_MonteCarlo_temp,type = "response")
        }
      }

      data_MonteCarlo_temp[,treat] <- a1
      y_hat <- predict(modelY,data_MonteCarlo_temp,type = "response")
      phi[p] <- mean(y_hat)
    }
    res <- t(as.matrix(c(phi[med.num+1]-phi0, phi[1:med.num]-phi[med.num+1], phi1-sum(phi[1:med.num])+(med.num-1)*phi[med.num+1],
                         phi1-phi0)))
    colnames( res ) <- c(paste("PSE",0:med.num,sep=""),"FMI","TE")
    return(FMI=  res )
  }

  set.seed(seed)
  res0 <- FMI_main(data=data,outcome = outcome,treat = treat,
                  mediators=mediators, baseline.conf = baseline.conf,
                  med.model=med.model,
                  outcome.model=outcome.model,
                  a0 = a0,a1 = a1, MonteCarlo.time = MonteCarlo.time)

  res <- matrix(NA,nboost,length(mediators)+3)
  for(b in 1:nboost){
    data.b <- data[sample(1:nrow(data),replace = T),]
    res[b,] <- FMI_main(data=data.b,outcome = outcome,treat = treat,
                    mediators=mediators, baseline.conf = baseline.conf,
                    med.model=med.model,
                    outcome.model=outcome.model,
                    a0 = a0,a1 = a1, MonteCarlo.time = MonteCarlo.time)
  }
  output <- matrix(NA,length(mediators)+3,4)
  output[,1] <-  res0
  output[,2] <-  apply(res, 2, sd)
  output[,3] <- output[,1]-1.96*output[,2]
  output[,4] <- output[,1]+1.96*output[,2]

  rownames( output ) <- c(paste("PSE",0:length(mediators),sep=""),"FMI","TE")
  colnames( output ) <- c("Estimate","SD","95%CI (lower bound)","95%CI (upper bound)")
  return(FMI=  output )
}
