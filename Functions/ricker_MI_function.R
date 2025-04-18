########################################################################
# This user defined function will fill in missing data with multiple imputation using Amelia
# and then fit the Ricker Model using Poisson or Negative Binomial error distribution


# Load packages
#library(here)
library(stats)
library(Amelia)
library(R.utils)
library(MASS)


#' Title Fit Ricker Model to population count data with Multiple Imputation using Amelia
#' 
#' This function
#'
#' @param y vector of population time series data with NAs in any positions 
#' @param imputationsnum number of imputations to average together, more will take longer but may have greater accuracy 
#' @param fam either fit the models using poisson method ("poisson") or negative binomial method ("neg_binom")
#' @param method Each missing data point will be filled in twice by Amelia, once as y_t and once as y_t-1, should Amelia only predict y_t ("forward"), only predict y_t-1 ("backward"), use both and allow discrepancies ("dual", default), or use both and average out discrepancies ("averaging")
#' @param p2samelia an integer value to pass to amelia taking either 0 for no screen output, 1 for normal screen printing of iteration numbers, and 2 for detailed screen output
#' @param ameliatimeout a number in seconds indicating when to cut off Amelia's fitting algorithm and indicate convergence issues, default is 1 minute
#' @param pro_conf what type of projection confidence to provide, defaults to "none", "sim" returns 1000 simulated estimates, "boot" returns 1000 bootstrapped estimates
#'
#' @return List of estimates, standard errors, and confidence intervals, or NA if an error occurred
#'
#' @examples
#' 
#' y <- readRDS("data/missingDatasets/pois_sim_randMiss_B.rds")[[1]]$y[[10]]
#' fit_ricker_MI(y,ameliatimeout=10)
#' 
fit_ricker_MI<-function(y, imputationsnum=5, fam = "poisson", method="dual", p2samelia=1, ameliatimeout=60, pro_conf="none"){
  
  # Check for population extinction
  if(sum(y==0,na.rm=T)>0){
    warning("population extinction caused a divide by zero problem, returning NA")
    return(list(
      NA,
      cause = "population extinction"
    ))
  }
  
  # Check for NaN
  if(any(is.nan(y),na.rm=T)){
    warning("NaN found, recode missing data as NA, returning NA")
    return(list(
      NA,
      reason = "NaN found"
    ))
  }
  
  # Check for Inf
  if(any(is.infinite(y),na.rm=T)){
    warning("infinite population detected, recheck data returning NA")
    return(list(
      NA,
      reason = "population explosion"
    ))
  }
  
  # check that there actually are missing values 
  if(any(is.na(y)) == FALSE) {
    warning("cannot fit a model with Amelia, there are no missing values to impute!")
    return(list(
      NA, 
      reason = "no missing values"
    ))
  }
  
  # get length of time series
  n <- length(y)
  
  # Check for not enough information for amelia to fit imputation model (amelia returns collinearity error)
  if(length(which(is.na(y[2:n]-y[1:(n - 1)])))>=(n-2)){ # if there are no, or only 1 overlap of non-NAs for Amelia to use for MI
    warning("There are not enough non-missing sets y(t) and y(t-1)")
    return(list(
      NA,
      reason = "missingness limits"
    ))
  }
  
  # check to see if the missingness limit has been met (i.e. if there are fewer than 3 non-missing values??)
  # fail if trimmed time series is too small 
  if (sum(!is.na(y)) <=5) {
    warning("Time series with NAs dropped is too short! Model can't fit well")
    return(list(
      NA,
      reason = "ts too short"
    ))
  }
  # another missingness check... make sure that there are more than three consecutive pairs of real data
  if (sum(!is.na(y-lag(y))) <=3){
    warning("too few consecutive time points (three or less)")
  }
  
  # make data frame in prep for multiple imputation
  simmissingdf=cbind.data.frame(1:(n-1),y[2:n],y[1:(n - 1)])
  names(simmissingdf)<-c("time","yt","yt1")
  
  
  # bounds on amelia guesses- population cannot be lower than 1 (negative numbers are a problem as are magical recoveries from extinctions)
  bound1=matrix(c(2,1,999999,3,1,999999),byrow=T,ncol=3)
  
  
  # Do MI with Amelia, with a time limit so that it can timeout 
  tryCatch(
    expr = {
      withTimeout(expr={
        # in the "dual" method we use amelia to get rid of NAs in both y(t) and y(t-1), without regard to values that *would* be the same, this method is the fastest
        if(method=="dual"){
          amelia1sim<-amelia(simmissingdf, m=imputationsnum, ts="time", bounds=bound1, p2s=p2samelia, leads="yt")
          for(i in 1:imputationsnum){
            while(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0|length(which(is.na(amelia1sim$imputations[[i]]$yt1)))>0){
              # get where the NAs are
              navec<-which(is.na(amelia1sim$imputations[[i]]$yt))
              # replace some NAs with values that should be the same
              amelia1sim$imputations[[i]]$yt[navec]<-amelia1sim$imputations[[i]]$yt1[navec+1]
              amelia1sim$imputations[[i]]$yt1[navec]<-amelia1sim$imputations[[i]]$yt[navec-1]
              # do another amelia round to fill in any more NAs if possible
              if(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
                am2<-amelia(amelia1sim$imputations[[i]],m=1,ts="time", bounds=bound1, p2s=p2samelia, leads="yt")
                amelia1sim$imputations[[i]]<-am2$imputations$imp1
              }
            }
          }
        }
        
        # in the "forward" method we use amelia to get rid of NAs in y(t), and fill in y(t-1) based on which population values should match
        if(method=="forward"){
          amelia1sim<-amelia(simmissingdf, m=imputationsnum, ts="time", bounds=bound1,p2s=p2samelia, leads="yt")
          for(i in 1:imputationsnum){
            #correct it to only predict forward, that is only predict yt, and then fill the yt1
            fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
            amelia1sim$imputations[[i]][fill1+1,3]=amelia1sim$imputations[[i]][fill1,2]
            while(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
              # get where the NAs are
              navec<-which(is.na(amelia1sim$imputations[[i]]$yt))
              # replace some NAs with values that should be the same
              #amelia1sim$imputations[[i]]$yt[navec]<-amelia1sim$imputations[[i]]$yt1[navec+1]
              amelia1sim$imputations[[i]]$yt1[navec]<-amelia1sim$imputations[[i]]$yt[navec-1]
              # do another amelia round to fill in any more NAs if possible
              if(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
                am2<-amelia(amelia1sim$imputations[[i]],m=1,ts="time", bounds=bound1,p2s=p2samelia, leads="yt")
                amelia1sim$imputations[[i]]<-am2$imputations$imp1
              }
              #correct it to only predict forward, that is only predict yt, and then fill the yt1
              fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
              amelia1sim$imputations[[i]][fill1+1,3]=amelia1sim$imputations[[i]][fill1,2]
              
            }
          }
        }
        
        # in the "backward" method we use amelia to get rid of NAs in y(t-1), and fill in y(t) based on which population values should match
        if(method=="backward"){
          amelia1sim<-amelia(simmissingdf, m=imputationsnum, ts="time", bounds=bound1,p2s=p2samelia, leads="yt")
          for(i in 1:imputationsnum){
            #correct it to only predict backward, that is only predict yt1, and then fill the yt
            fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
            amelia1sim$imputations[[i]][fill1,2]=amelia1sim$imputations[[i]][fill1+1,3]
            while(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
              # get where the NAs are
              navec<-which(is.na(amelia1sim$imputations[[i]]$yt))
              # replace some NAs with values that should be the same
              amelia1sim$imputations[[i]]$yt[navec]<-amelia1sim$imputations[[i]]$yt1[navec+1]
              #amelia1sim$imputations[[i]]$yt1[navec]<-amelia1sim$imputations[[i]]$yt[navec-1]
              # do another amelia round to fill in any more NAs if possible
              if(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
                am2<-amelia(amelia1sim$imputations[[i]],m=1,ts="time", bounds=bound1,p2s=p2samelia, leads="yt")
                amelia1sim$imputations[[i]]<-am2$imputations$imp1
              }
              #correct it to only predict backward, that is only predict yt1, and then fill the yt
              fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
              amelia1sim$imputations[[i]][fill1,2]=amelia1sim$imputations[[i]][fill1+1,3]
              
            }
          }
        }
        
        # in the "averaging" method we use amelia to get rid of NAs in y(t) and y(t-1), then average values that have been predicted twice
        if(method=="averaging"){
          amelia1sim<-amelia(simmissingdf, m=imputationsnum, ts="time", bounds=bound1,p2s=p2samelia, leads="yt")
          for(i in 1:imputationsnum){
            
            while(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
              # get where the NAs are
              navec<-which(is.na(amelia1sim$imputations[[i]]$yt))
              # replace some NAs with values that should be the same
              amelia1sim$imputations[[i]]$yt[navec]<-amelia1sim$imputations[[i]]$yt1[navec+1]
              amelia1sim$imputations[[i]]$yt1[navec]<-amelia1sim$imputations[[i]]$yt[navec-1]
              # do another amelia round to fill in any more NAs if possible
              if(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
                am2<-amelia(amelia1sim$imputations[[i]],m=1,ts="time", bounds=bound1,p2s=p2samelia, leads="yt")
                amelia1sim$imputations[[i]]<-am2$imputations$imp1
              }
            }
            #correct it to average the predicted yt and yt1 that don't match
            fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
            fill2=(amelia1sim$imputations[[i]][fill1,2]+amelia1sim$imputations[[i]][fill1+1,3])/2
            amelia1sim$imputations[[i]][fill1,2]=fill2
            amelia1sim$imputations[[i]][fill1+1,3]=fill2
          }
        }
        
      },timeout = ameliatimeout)
    },
    TimeoutException=function(msg){
      warning("Amelia has timed out, likely due to high missingness")
      return(list(
        NA,
        reason = "Amelia time out"
      ))
    }, error=function(e){
      warning("Amelia unable to fit the model, likely due to too little data")
      return(list(
        NA,
        reason="Amelia fitting error"
      ))
    }
  )
  
  if(!exists("amelia1sim")){
    warning("Amelia was unable to fit for reason other than timeout")
    return(list(
      NA,
      reason = "Amelia fitting error"
    ))
  } 
  
    if(amelia1sim$code!=1){
      cat("we should be returning NA, code is not 1")
      warning("Amelia was unable to fit for reason other than timeout")
      return(list(
      NA,
      reason = paste("Amelia internal fitting error, code",amelia1sim$code)
    ))
  } 
  
  
  if(any(is.na(amelia1sim))){
    warning("Amelia has timed out, likely due to exceptionally high missingness")
    return(list(
      NA,
      reason = "Amelia time out"
    ))
  } 
  
  if(any(is.na(amelia1sim$imputations))){
    warning("Amelia has timed out, likely due to exceptionally high missingness")
    return(list(
      NA,
      reason = "Amelia fitting error"
    ))
  } 
  
  
  fit=list()
  # fit model over all imputations
  for(i in 1:imputationsnum){
    
    # compile into sliced dataframe in preparation for 
    dat =data.frame(
      yt = amelia1sim$imputations[[i]][2:(n-1),2],
      ytm1 = amelia1sim$imputations[[i]][1:(n - 2),2]
    )
    
    
    # fit ricker model with poisson
    if(fam == "poisson"){
      fit[[i]]<-tryCatch({
        glm(
          round(yt) ~ round(ytm1), data = dat, 
          family = poisson,
          offset = log(ytm1)
        )
      },error=function(cond){
        message(paste("we have had an error in the model fitting X2"))
        return(list(NA, reason="model fitting error"))
      }
      )
      
    }
    
    # or fit ricker model with negbinom
    if(fam == "neg_binom"){
      fit[[i]] <- tryCatch({
        MASS::glm.nb(
          round(yt) ~ round(ytm1) + offset(log(round(ytm1))), 
          data = dat
        )
      },error=function(cond){
        message(paste("we have had an error in the model fitting"))
        return(list(NA, reason="model fitting error"))
      }) 
    }
    
  }
  
  # Check that the model fit ran ok
  if(any(is.na(fit))){
    warning("There has been a model fitting error")
    return(list(NA,
                reason="Amelia model fit error"
                  ))
  }
  
  
  sapply(fit, profile)
  # Averages models into 1 result, simplifies, renames
  estims1=sapply(fit,coef,simplify=T)
  estims=rowMeans(estims1)
  estims=estims * c(1, -1)
  names(estims) <- c("r", "alpha")
  
  cis1=sapply(fit,confint,simplify=T)
  cis=rowMeans(cis1)
  l1=as.double(cis[c(1,4)] * c(1, -1))
  names(l1) <- c("r", "alpha")
  u1=as.double(cis[c(3,2)] * c(1, -1))
  names(u1) <- c("r", "alpha")
  
  ses1=sapply(fit,function(X) sqrt(diag(vcov(X))),simplify=T)
  ses=rowMeans(ses1)
  names(ses) <- c("r", "alpha")
  
  # no need for 1000 betas for projection CI
  if(pro_conf=="none"){
    # return as a list
    return(list(
      estim = estims,
      se = ses,
      lower = l1,
      upper = u1
    ))
  }
  
  # simulated 1000 betas for projection CI
  if(pro_conf=="sim"){
    sim_per_fit=1000/length(fit)

    results <- lapply(fit, function(model) {
      betas <- coef(model)     # Get coefficients
      vcov_mat <- vcov(model)  # Get variance-covariance matrix
      mvrnorm(sim_per_fit, betas, vcov_mat)  # Generate random samples
    })

    CI_results=do.call(rbind, results)
    
    return(list(
      estim = estims,
      se = ses,
      lower = l1,
      upper = u1,
      CI_results=CI_results
    ))
    
  }
  
  # bootstrapped 1000 betas for projection CI
  if(pro_conf=="boot"){
    sim_per_fit=1000/length(fit)
    
    boot = function(x, model){
      data = model.frame(model)
      data_sample = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      names(data_sample)=c("yt","ytm1","offset")
      coef(update(model, data = data_sample))
    }
    
    CI_results <- do.call(rbind, lapply(fit, function(single_fit) {
      do.call(rbind, lapply(1:sim_per_fit, boot, single_fit))
    }))
    
    return(list(
      estim = estims,
      se = ses,
      lower = l1,
      upper = u1,
      CI_results=CI_results
    ))
  }
  
}


