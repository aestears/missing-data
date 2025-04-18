
# Load packages -----------------------------------------------------------

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(plyr)
library(shinystan)
library(faux)
library(bayesplot)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(data.table)
library(Amelia)

# Simulate data -----------------------------------------------------------

N<-365 #length of data
z<-numeric(N+1) 
t<-numeric(N+1) 
time<-seq(as.POSIXct("2020/1/1"),as.POSIXct("2020/12/30"), by="day") #for light



##light
# From Yard et al. (1995) Ecological Modelling.  Remember your trig?  
# calculate light as umol photon m-2 s-1.
# Arguments are:  
# time = a date and time input (posixct object)
# lat = latitude of field site
# longobs = longitude of field site
# longstd = standard longitude of the field site (NOTE: watch daylight savings time!!!). For PST, longstd is be 120 degrees. But during PDT it is 105 degrees. MST is 105 deg. MDT is 90. 


# convert degrees to radians
radi<-function(degrees){(degrees*pi/180)}

# function to estimate light
lightest<- function (time, lat, longobs, longstd) {
  jday<-yday(time)
  E<- 9.87*sin(radi((720*(jday-81))/365)) - 7.53*cos(radi((360*(jday-81))/365)) - 1.5*sin(radi((360*(jday-81))/365))
  LST<-as.numeric(time-trunc(time))
  ST<-LST+(3.989/1440)*(longstd-longobs)+E/1440
  solardel<- 23.439*sin(radi(360*((283+jday)/365)))
  hourangle<-(0.5-ST)*360
  theta<- acos( sin(radi(solardel)) * sin(radi(lat)) + cos(radi(solardel)) * cos(radi(lat)) * cos(radi(hourangle)) )
  suncos<-ifelse(cos(theta)<0, 0, cos(theta))
  GI<- suncos*2326
  GI	
  
}
light<-lightest(time, 47.8762, -114.03, 105) #Flathead Bio Station just for fun
light.rel<-light/max(light)


#Load fake Q
setwd("~/GitHub/missing-data")
discharge<-read_csv("fake_turb.csv")
q.rel<-discharge$turb/max(discharge$turb)


##GPP based on time-series model with known parameters
sdp <- 0.1
sdo<-0.01
phi <-0.8
b1<-0.6
b2<--0.8
b0<-2

# center at expected value:
EV<- (b0*(1-phi)+(b1*light.rel[1])+(b2*q.rel[1]))/(1-phi) #Expected value ts[1] varies with variable shoice
z[1]<-EV
## Set the seed, so we can reproduce the results
# set.seed(5)
## For-loop that simulates the state through time, using i instead of t,
# In this process equation, I have centered the data to help remove the equifinality
# between b0 and phi. See: 
# https://stats.stackexchange.com/questions/348430/what-is-the-reason-for-not-including-an-intercept-term-in-ar-and-arma-models
# with no intercept: y_t = phi * y_t-1 + b1*light_t + w_t
# center y by subtracting the expected intercept (b0):
#                   y_t - b0 = phi*(y_t-1 - b0) + b1 * light_t + w_t
# rearrange:        y_t = b0 + phi * y_t-1 - phi*b0 + b1*light_t + w_t
#                   y_t = b0*(1 - phi) + phi * y_t-1 + b1*light_t + w_t
# now with the centered data, the intercept is defined in terms of phi 
for (i in 1:365){
  z[i+1] = b0 * (1 - phi) + z[i] * phi + light.rel[i] * b1 + q.rel[i]*b2 +rnorm(1, 0, sdp)
}

# Estimate observed data
# set.seed(4)
sd.o <-rnorm(N, 0, sdo)
y <-z[2:(N+1)] + sd.o
y[N+1]<-NA

##Bind simulated data and time
y_full<-as.data.frame(cbind(z,y,seq(0:N+1)))

ggplot(data=y_full, aes(x=V3, y=z))+
  geom_point(aes(color="state"), size=2)+
  geom_point(aes(y=y, x=V3,color="observed"), size=2)+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  theme(legend.position="top")+
  ylab("Data")+
  xlab("Time (days)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_color_manual(name=c("observed","state"), values = c("red","black"))


# Force some missing data-BY DAY ------------------------------------------
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created

##Set seed for replicability
set.seed(50)

##Load data from SDW characterizing lengths
obs<-readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/SDW_obs_length.RDS")
miss<-readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/SDW_miss_length.RDS")

##Replicate simulated data
n.datasets<-41
y_miss<-list(y)
y_miss<-rep(y_miss[1], times=n.datasets)


##set up loop to create missing datasets from y_miss. 
#The loop pulls random lengths of observed and missing data and makes missing data zero. 
#It switches between pulling a missing length and replaces those values with zero,
#then it picks an observed length and jumps ahead (keeping the data points in the jump).
#The loop for each individual dataset ends when it exceeds the length N
#The end result are N=n.datasets with varying missing data gaps and percents of missing data

#Create objects for the loop
sample.miss<-list(NA)
sample.obs<-list(NA)
sample.miss<-rep(sample.miss, times=n.datasets)
sample.obs<-rep(sample.obs, times=n.datasets)

#For loop for N=n.datasets. The end result are N=n.datasets with 0 for missing data
for (i in 2:n.datasets){
#samples 200 times from SDW lengths with replacement
sample.miss[[i]]<-sample(miss,200,replace=TRUE)
sample.obs[[i]]<-sample(obs,200,replace=TRUE)

#Initiates new dataset. It will alwasys start with some length of observed data. 
new<-sample.obs[[i]][1]

#For loop for each individual dataset
for(j in 1:365){
  
  remove<-seq(from=new+1, to=new+sample.miss[[i]][j], by=1)
  y_miss[[i]][remove]<-0
  new<-new+sample.miss[[i]][j]+sample.obs[[i]][j+1]
  if (new>=N){
    break
  }
}
}

#Assign NAs to the missing data in each list
for(i in 1:n.datasets){
  z<-y_miss[[i]]
  z[which(z==0)]<-NA
  y_miss[[i]]<-z
}

#Cuts each down to size, as they will often end with missing data
for (i in 1:n.datasets){
  y_miss[[i]]<-y_miss[[i]][1:365]
}


#Check to make sure
summary(y_miss)


#Calculate the average length of data gaps in each dataset
## Calculate length of data gaps
# Prepare vectors for loop
miss_vec<-list(NA)
miss_vec<-rep(miss_vec, times=n.datasets)
length_vec<-list(NA)
length_vec<-rep(length_vec, times=n.datasets)

# Create dummy variable for present or not
dummy<-list(NA)
dummy<-rep(dummy, times=n.datasets)

for (i in 1:n.datasets){
  dummy[[i]]<-ifelse(is.na(y_miss[[i]]),1,0)
}

#loop over data and return the length of missing data gaps
for (i in 1:n.datasets){
  # First data point will always be observed
  miss_vec[[i]][1]<-dummy[[i]][1]
  length_vec[[i]][1]<-1
  
  for (j in 1:N){
    miss_vec[[i]][j]<-ifelse(dummy[[i]][j]==1,1+miss_vec[[i]][j-1], 0) #If data is missing add to yesterdays value by 1. If data observed then 0.
    length_vec[[i]][[j]]<-ifelse(miss_vec[[i]][j]==0,miss_vec[[i]][j-1], NA) #Save length of gap on the day data is observed (0) after a gap
    # length_miss[i]<-ifelse(length_miss[i]==0,NA,length_miss[i]) #Replace days with observed data with NA
  }
}

##Calculate mean of missing data gaps
missing<-lapply(length_vec,na.omit)
missing<-lapply(missing, function(missing) {missing[missing!=0]})
mean.missing.gap<-unlist(lapply(missing,mean))
mean.missing.gap

##first dataset has no gaps (is NaN)
mean.missing.gap[[1]]<-0
hist(mean.missing.gap)
#Calculate the percent of missing data
miss.num<-sapply(y_miss, function(y) sum(length(which(is.na(y)))))
prop.miss<-round(miss.num/length(y)*100)
prop.miss
hist(prop.miss)

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
#FOR SINGLE VECTOR
#y_index_mis <- which(is.na(y_miss)) # Identify the rows for each variable in airquality with NA
#y_index_obs <- which(!is.na(y_miss)) # Identify the rows for each variable in airquality with observed data
#y_nMiss <- length(y_index_mis)# How many NAs?
#y_nObs <- length(y_index_obs) # How many NAs?

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
#FOR LIST OF VECTORS
y_index_mis <- lapply(y_miss,function(var){which(is.na(var))}) # Identify the rows for each variable in airquality with NA
y_index_obs <- lapply(y_miss,function(var){which(!is.na(var))}) # Identify the rows for each variable in airquality with observed data
y_nMiss <- lapply(y_index_mis,function(var){length(var)}) # How many NAs?
y_nObs <- lapply(y_index_obs,function(var){length(var)}) # How many NAs?


##Replace NAs with arbitrary number to make Stan happy
for(i in 1:n.datasets){
  r<-y_miss[[i]]
  r[y_index_mis[[i]]]<--100 #arbitrary number
  y_miss[[i]]<-r
} 

summary(y_miss)



# Model code for STAN -----------------------------------------------------


#Stan model

#see https://michaeldewittjr.com/resources/stan_nowcast.html#

sink("ss_light_centered_miss.stan")

cat("

  
/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] y; // Observations
    vector[N] light;
    vector[N] q;
    real z0;
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss];
    vector[N] miss_vec; //index or location of missing values within the dataset
    }
  
 
/*----------------------- Parameters --------------------------*/
  
  parameters {
    vector[N] z; //latent state variable
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real b0; // intercept
    real b1; // light coefficient 
    real b2; //Q coefficient
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
      }
 
 
  /*----------------------- Model --------------------------*/
  model {
   
  //merge imputed and observed data
    // Prior distributions
    sdo ~ normal(0.01, 0.001);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    b1 ~ normal(0,5);
    b2 ~ normal(0,5);
    
    
   // Distribution for the zero and first state
    
   // Distribution for the first state
  z[1] ~ normal(z0, sdp);
   y[1]~ normal(z[1],sdo);
   
    // Process model
   
    for(t in 2:N){
      z[t] ~ normal(b0*(1-phi)+ z[t-1]*phi+light[t]*b1+ q[t]*b2, sdp);
    if(miss_vec[t]==0){
        target += normal_lpdf(y[t]| z[t], sdo); // observation model with fixed observation error
    }
    }
    
 }
   
      "
    
    ,fill=TRUE)
sink()
#closeAllConnections()





# DA-model (with loop) -----------------------


fit.miss <- vector("list",n.datasets)
model<-"ss_light_centered_miss.stan"
model<-stan_model(model)
for(i in 2){
##Load data
  data <- list(   N = length(y_miss[[i]]),
                  y_nMiss = unlist(y_nMiss[[i]]),
                  y_index_mis = unlist(y_index_mis[[i]]),
                  y= unlist(y_miss[[i]]),
                  light=light.rel,
                  q=q.rel[1:365],
                  z0=y_miss[[i]][1],
                  miss_vec=miss_vec[[i]])
  fit.miss[[i]]<-rstan::sampling(object=model, data = data,  iter = 3000, chains = 4, control=list(max_treedepth = 15))
  print(i)
  rstan::check_hmc_diagnostics(fit.miss[[i]])
}


# DA-pull and clean parameters ---------------------------------------------

traceplot(fit.miss[[1]], pars=c("sdo","sdp", "phi", "b0", "b1","b2"))
##Pull param estimates into list
fit_summary_pars_bayes <- vector("list",41)
for (i in 1:41){
   fit_summary_pars_bayes[[i]]<-(summary(fit.miss[[i]], pars=c("sdp","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary)
   }
fit_summary_pars_bayes[[1]]

##Unlist,cleanup, and add factors
fit_summary_bayes<-as.data.frame(do.call(rbind.data.frame, fit_summary_pars_bayes)) #Unlist
fit_summary_bayes$param<-rep(c("sdp","phi", "b1","b0"), times=41 )#add parameter factor
fit_summary_bayes$prop.missing<-rep(prop.miss[1:41], times=rep(4, times=41)) #add prop of missing data
fit_summary_bayes$gap.missing<-rep(mean.missing.gap[1:41], times=rep(4, times=41)) #add prop of missing data
row.names(fit_summary_bayes)<-1:(41*4) #remove row names
fit_summary_bayes$param<-as.factor(fit_summary_bayes$param) #make factor
#fit_summary_bayes$prop.missing<-as.factor(fit_summary_bayes$prop.missing) #make factor
colnames(fit_summary_bayes)<-c("mean", "se-mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing", "gap.missing")
summary(fit_summary_bayes) #check it looks good
head(fit_summary_bayes)

fit_summary_bayes<-fit_summary_bayes[order(fit_summary_bayes$prop.missing),]
known.data<-c(sdp,phi,b1,b0)
known<-as.data.frame(known.data)

known.param<-c("sdp","phi", "b1","b0")
#known.missing<-c(5, 5, 5,5)
known<-as.data.frame(cbind(known.data, known.param))
known$known.data<-as.numeric(as.character(known$known.data))

ggplot(data=fit_summary_bayes, aes(x=mean, y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3)+
  theme_classic()+
  theme(legend.position="top")+
  ylab("Percent of Missing Data")+
  xlab("Parameter Estimate")+
  scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
  #scale_x_continuous(limits=c(0,1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept = known.data,color=c("black", "darkgray", "green", "blue"))



# DA-Write or read RDS file --------------------------------------------------

write_rds(fit_summary_bayes,"C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/summary_sim_var_bayes_sdp_01_phi_8_b0_1_b1_1.RDS")

fit_summary_bayes <- readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/summary_sim_var_bayes_sdp_1_phi_8_b0_1_b1_1.RDS")



# DA-calculate RWCI and bias ----------------------------------------------


##Calculate denominator of rwci for easier calculation
fit_summary_bayes$rwci.den<-rep(fit_summary_bayes$high[1:4]-fit_summary_bayes$min[1:4], times=41)
fit_summary_bayes$known.param<-rep(c(0.1, 0.8, 0.1, 0.1), times=41)

##Calculate bias (difference from known) and rwci (relative width of the credible interval)

diff<-ddply(fit_summary_bayes, c("prop.missing", "param"),summarize,
            bias=mean-known.param, 
            rwci= (high-min)/rwci.den
)

fit_summary_bayes$bias<-diff$bias
fit_summary_bayes$rwci<-diff$rwci
known.param<-as.data.frame(cbind(as.numeric(fit_summary_bayes$known.param[1:4]), as.character(fit_summary_bayes$param[1:4])))
colnames(known.param)<-c("known","param")
known.param$known<-as.numeric(known.param$known)

ggplot(data=fit_summary_bayes, aes(x=mean, y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3)+
  theme_classic()+
  theme(legend.position="right")+
  xlab("Parameter value")+
  ylab("Percent of Missing Data")+
  scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
  #scale_x_continuous(limits=c(-1,1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept =known.param$known,color=c("black","darkgray", "blue", "green" ))

ggplot(data=fit_summary_bayes, aes(x=log(rwci), y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3)+
  theme_classic()+
  theme(legend.position="right")+
  ylab("Percent of Missing Data")+
  xlab("log Relative width of 95% credible interval log")+
  scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
  #scale_x_continuous(limits=c(-.05,1.1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept = 0,color=c("black"))

ggplot(data=fit_summary_bayes, aes(x=bias, y=gap.missing ))+
  geom_point(aes( color=param, group=param),size=3)+
  theme_classic()+
  theme(legend.position="right")+
  xlab("Parameter value")+
  ylab("Percent of Missing Data")+
  scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
  #scale_x_continuous(limits=c(-1,1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept =known.param$known,color=c("black","darkgray", "blue", "green" ))




# DA-model check ----------------------------------------------------------


##Print and extract
fit_extract<-rstan::extract(fit.miss[[2]])
print(fit.miss[[2]], pars=c( "sdp","phi", "b1", "sdo","b0"))

pairs(fit, pars=c("sigma_proc","phi", "b1", "sigma_obs"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit.miss[[16]])

##Traceplots

traceplot(fit.miss[[1]], pars=c("phi", "b0","sdp", "b1"))

##Posterior densities compared to known parameters
plot_sdo <- stan_dens(fit.miss[[2]], pars="sdo") + geom_vline(xintercept =sdo)
plot_sdp <- stan_dens(fit.miss[[2]], pars="sdp") + geom_vline(xintercept = sdp)
plot_phi <- stan_dens(fit.miss[[2]], pars="phi") + geom_vline(xintercept = phi)
plot_b2 <- stan_dens(fit.miss[[2]], pars="b2") + geom_vline(xintercept = b2)+xlab("Q beta")
plot_b1 <- stan_dens(fit.miss[[2]], pars="b1") + geom_vline(xintercept = b1)+xlab("Light beta")
plot_b0 <- stan_dens(fit.miss[[2]], pars="b0") + geom_vline(xintercept = b0)+xlab("intercept")
grid.arrange(plot_b1, plot_b2,plot_phi,plot_sdp,plot_sdo, plot_b0, nrow=2)

##Posterioir Predictive check and test statistic
yrep1 <- fit_extract$y_rep
samp100 <- sample(nrow(yrep1), 100)
ppc_dens_overlay(y_miss[[6]], yrep1[samp100, ]) 
ppc_stat(y, yrep1[samp100, ])

##Shineystan for everything else
launch_shinystan(fit)

##Create object with estimated missing data

fit_summary<-summary(fit.miss[[2]], probs=c(0.05,.5,.95))$summary
my_data <- as_tibble(fit_summary)
y_est_data<-my_data[1:365,]


##Merge observed and estimated

y_combined<-cbind(y_full$z[2:366],y_full$y[1:365],y_full$V3[1:365],y_est_data)

##rename column names for clarity

colnames(y_combined)<-c("obs.z", "obs.y", "time", "est.z", "se.z","sd.z", "low", "median","high", "n_eff", "rhat")
##rename column names for clarity

#colnames(y_combined)<-c("date", "observed", "estimated", "se_mean_est", "sd_est", "low", "median","high", "n_eff", "rhat")

##check that column names are correct

#head(y_combined)

##Plot observed and estimated direct comparison
ggplot(data=y_combined, aes(x=obs.z, y=est.z))+
  geom_point( size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()+
  geom_errorbar(aes(ymin=low, ymax=high))+
  theme(legend.position="top")+
  ylab("Estimated GPP")+
  xlab("Observed GPP")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))


raw_data<-as.data.frame(cbind(z,time))

ggplot(data=y_combined, aes(x=date, y=estimated))+
  geom_point(data=raw_data, aes(x=date, y=z, color="gray"), size=3)+
  geom_errorbar(data=y_combined,aes(x=date,ymin=low, ymax=high))+ 
  geom_point(aes(color="red"), size=3)+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  scale_color_identity(guide = "legend", name=element_blank(), labels=c("Simulated data", "Imputed data"))+
  theme(legend.position="top")+
  ylab("Data")+
  xlab("Time")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

saveRDS(fit.miss, file = "full_sim_var_bayes_sdp_01_phi_8_b0_1_b1_1_part4.RDS")
saveRDS(fit_summary_bayes, file = "summary_sim_var_bayes_sdp_01_phi_8_b0_1_b1_1_allparts.RDS")
saveRDS(known, file = "known_sim_var_bayes_sdp_01_phi_8_b0_1_b1_1_part4.RDS")


ts<-as.Date(time,format="%m/%d/%Y")



# MI-create datasets ------------------------------------------------------


data.amelia<- vector("list",41)
for(i in 2:41){
y_miss1<-na_if(y_miss[[i]], -100)
y_miss1[1]<-z[1]
w<-as.data.frame(cbind(y_miss1,light,time))
z2<-amelia(w, m = 5, p2s=1, ts="time", lags="y_miss1")
summary(z2)
data.amelia[[i]]<-z2
}


#Assign NAs to the missing data in each list
for(i in 2:6){
  z1<-data.amelia[[i]]$imputations$imp1$y_miss1
  z2<-data.amelia[[i]]$imputations$imp2$y_miss1
  z3<-data.amelia[[i]]$imputations$imp3$y_miss1
  z4<-data.amelia[[i]]$imputations$imp4$y_miss1
  z5<-data.amelia[[i]]$imputations$imp5$y_miss1
  q1<-z1[y_index_mis[[i]]]
  q2<-z2[y_index_mis[[i]]]
  q3<-z3[y_index_mis[[i]]]
  q4<-z4[y_index_mis[[i]]]
  q5<-z5[y_index_mis[[i]]]
  q.ts<-time[y_index_mis[[i]]]
  q<-cbind(q1,q2,q3,q4,q5)
  q.mean<- apply(q[,-1], 1, mean)
  q.min<-apply(q[,-1], 1, min)
  q.max<-apply(q[,-1], 1, max)
  q.sd<- apply(q[,-1], 1, sd)
  x_obs_data<-as.data.frame(cbind(z,time))
  x_imp_data<-as_tibble(cbind(q.mean,q.sd,q.ts, q.min,q.max))
  x_imp_data<-as.data.frame(x_imp_data %>% slice(1:length(y_index_mis[[i]])))
  x_obs<-z[y_index_mis[[i]]]
  x_combined<-as.data.frame(cbind(x_imp_data$q.mean,x_imp_data$q.sd,x_obs,x_imp_data$q.min,x_imp_data$q.max))
  g<-ggplot(data=x_imp_data, aes(x=q.ts, y=q.mean))+
    geom_errorbar(data=x_imp_data,aes(x=q.ts,ymin=q.min, ymax=q.max))+
    geom_point(data=x_obs_data, aes(x=time, y=z, color="gray"), size=3)+
    geom_point(aes(color="red"), size=3)+
    theme_classic()+
    theme(axis.title.x=element_text(size=18,colour = "black"))+
    theme(axis.title.y=element_text(size=18,colour = "black"))+
    theme(axis.text.y=element_text(size=18,colour = "black"))+
    theme(axis.text.x=element_text(size=18,colour = "black"))+
    scale_color_identity(guide = "legend", name=element_blank(), labels=c("Simulated data", "Imputed data"))+
    theme(legend.position="top")+
    ylab("Data")+
    xlab("Time")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  s<-ggplot(data=x_combined, aes(x=x_obs, y=V1))+
    geom_abline(intercept=0, slope=1)+
    theme_classic()+
    geom_errorbar(aes(ymin=q.min, ymax=q.max))+
    geom_point( size=3)+
    theme(legend.position="top")+
    ylab("Imputed data (mean and range)")+
    xlab("Simulated data")+
    theme(axis.title.x=element_text(size=18,colour = "black"))+
    theme(axis.title.y=element_text(size=18,colour = "black"))+
    theme(axis.text.y=element_text(size=18,colour = "black"))+
    theme(axis.text.x=element_text(size=18,colour = "black"))
  compare.density(data.amelia[[i]], var = "y_miss1", xlab="Data", main="", lwd=2)
  plot(g)
  plot(s)
  
}
     


# MI-model (with loop) ----------------------------------------------------


fit.amelia <- vector("list",length(missing_n_week))
data.stan.amelia <- vector("list",(length(missing_n_week)))
list.1<- vector("list",(5))
for(i in 2:length(missing_n_week)){##two because the list holding the lists starts with a list of misc data that arent the estimates
  for(g in 1:5){
 list.1[[g]] <- data.amelia[[i]]$imputations[[g]]$y_miss1
         
}
  data.stan.amelia[[i]]<-list.1
  }

plot(data.stan.amelia[[6]][[2]],data.stan.amelia[[6]][[3]])

##Run Stan


  fit.stan.miss.amelia <- vector("list",length(missing_n_week))
  fit.stan.miss.imp<- vector("list",5)
  model<-"toy2p.stan"
  model<-stan_model(model)
  for(i in 4:6){
    for(g in 1:5){
    ##Load data
    data <- list(  N = length(y_miss[[i]]),
                   y_nMiss = unlist(y_nMiss[[1]]),
                   y_index_mis = unlist(y_index_mis[[1]]),
                   y_miss= unlist(data.stan.amelia[[i]][[g]]),
                   light=light.l,
                   z0=y_miss[[i]][1])
    
    ##Run Stan
    fit.stan.miss.imp[[g]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
  }
    fit.stan.miss.amelia[[i]]<-fit.stan.miss.imp
  }


# MI-pull and clean paramters ---------------------------------------------

 
  ##Pull param estimates into list
  fit_summary_pars_amelia <- vector("list",6)
  list.2<- vector("list",6)
  for (i in 2:6){
    for (g in 1:5){
      list.2[[g]] <-summary(fit.stan.miss.amelia[[i]][[g]], pars=c("sdp","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary
  }
  fit_summary_pars_amelia[[i]]<-list.2
  }
  fit_summary_pars_amelia[[2]]
  
  ##Unlist,cleanup, and add factors
  fit_summary<- vector("list",6)
  for(i in 2:6){
    fit_summary[[i]]<-as.data.frame(do.call(rbind,fit_summary_pars_amelia[[i]]))
    fit_summary[[i]]$param<-rep(c("sdp","phi", "b1","b0"), times=length(5) )#add parameter factor
    fit_summary[[i]]$prop.missing<-rep(prop.miss[i], each=4) #add prop of missing data
    row.names(fit_summary[[i]])<-1:20 #remove row names
    fit_summary[[i]]$param<-as.factor(fit_summary[[i]]$param) #make factor
    fit_summary[[i]]$prop.missing<-as.factor(fit_summary[[i]]$prop.missing) #make factor
    colnames(fit_summary[[i]])<-c("mean", "se-mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing")
  }
  fit_summary<-as.data.frame(do.call(rbind,fit_summary))
  ##Unlist,cleanup, and add factors
  fit_summary_bayes<-as.data.frame(do.call(rbind.data.frame, fit_summary_pars_bayes)) #Unlist
  fit_summary_bayes$param<-rep(c("sdp","phi", "b1","b0"), times=1 )#add parameter factor
  fit_summary_bayes$prop.missing<-rep(0, times=c(4)) #add prop of missing data
  row.names(fit_summary_bayes)<-1:4 #remove row names
  fit_summary_bayes$param<-as.factor(fit_summary_bayes$param) #make factor
  fit_summary_bayes$prop.missing<-as.factor(fit_summary_bayes$prop.missing) #make factor
  colnames(fit_summary_bayes)<-c("mean", "se-mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing")
  summary(fit_summary_bayes) #check it looks good
  head(fit_summary_bayes)
  
  fit_summary<-rbind(fit_summary, fit_summary_bayes[1:4,])
  fit_summary$prop.missing<- factor(fit_summary$prop.missing, levels=c(0,8,13,25,35,59))
  summary(fit_summary) #check it looks good
  head(fit_summary)
  
  
  summ<-ddply(fit_summary, c("param", "prop.missing"), summarize,
              avg    = mean(mean),
              ci.min = mean(min),
              ci.max = mean(high))
  
  #,vw= sum(se.mean.sq)/5)
  
  fit_summary<-left_join(fit_summary, summ, by=c("param", "prop.missing"))
  
  
  known.data<-c(sdp,phi,b1,b0)
  known<-as.data.frame(known.data)
  
  known.param<-c("sdp","phi", "b1","b0")
  #known.missing<-c(5, 5, 5,5)
  known<-as.data.frame(cbind(known.data, known.param))
  known$known.data<-as.numeric(as.character(known$known.data))
  
  ggplot(data=fit_summary, aes(x=avg, y=prop.missing ))+
    geom_point(aes( color=param, group=param),size=3,position=position_dodge(0.5))+
    #geom_point(data=fit_summary[[3]],aes(x=mean, y=prop.missing, color=param, group=param),size=3,position=position_dodge(0.5))+
    #geom_point(data=fit_summary[[4]],aes(x=mean, y=prop.missing, color=param, group=param),size=3,position=position_dodge(0.5))+
    theme_classic()+
    geom_errorbar(aes(xmin=ci.min, xmax=ci.max,color=param, group=param), width=0.2,position=position_dodge(0.5))+
    theme(legend.position="top")+
    ylab("Percent of Missing Data")+
    xlab("Parameter Estimate")+
    scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
    scale_x_continuous(limits=c(0,0.85))+
    theme(axis.title.x=element_text(size=18,colour = "black"))+
    theme(axis.title.y=element_text(size=18,colour = "black"))+
    theme(axis.text.y=element_text(size=18,colour = "black"))+
    theme(axis.text.x=element_text(size=18,colour = "black"))+
    geom_vline(xintercept = known.data,color=c("black", "darkgray", "green", "blue"))


# MI-write or read RDS files ----------------------------------------------

  saveRDS(fit.stan.miss.amelia, file = "full_sim_week_amelia_sdp_01_phi_8_b0_1_b1_1.RDS")
  saveRDS(fit_summary, file = "summary_sim_week_amelia_sdp_01_phi_8_b0_1_b1_1.RDS")
  saveRDS(known, file = "known_sim_week_amelia_sdp_01_phi_8_b0_1_b1_1.RDS")
  
  
  
  ####OLD
  
  sink("mcelreath_miss_ss.stan")
  
  cat("
 

/*----------------------- Functions --------------------------*/  
  functions{
    vector merge_missing( int[] y_index_mis, vector z, vector y_imp) {
    int N = dims(z)[1];
    int N_miss = dims(y_imp)[1];
    vector[N] merged;
    merged = z;
    for (i in 1:N_miss)
        merged[y_index_mis[i] ] =y_imp[i];
    return merged;    
    }
  }
  
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int N; // Length of state and observation time series
    vector[N] y_miss; // Observations
    real z0; // Initial state value
    vector[N] light;
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss];
    vector[N] miss_vec// index or location of missing values within the dataset
  }
  
/*----------------------- Parameters --------------------------*/
  /* Parameter block: defines the variables that will be sampled */
  parameters {
    vector<lower = 0>[y_nMiss] y_imp;// Missing data
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real b0;
    real b1;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    vector[N] z; // State time series
  }
  

  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
    //merge imputed and observed data
    vector[N] B_merge;
    B_merge= merge_missing(y_index_mis, to_vector(z), y_imp);
   
    // Prior distributions
    sdo ~ normal(0, 1);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    b1 ~ normal(0,5);
    
    // Distribution for the first state
   B_merge[1] ~ normal(z0, sdp);
   y_miss[1]~ normal(B_merge[1],sdo);
   
    // Distributions for all other states
    for(t in 2:N){
       B_merge[t] ~ normal(b0+ B_merge[t-1]*phi+light[t]*b1, sdp);
    if(i%miss_vec==0){
      count = count + 1;
      target += normal_lpdf(y_miss[count]| B_merge[i], sdo); // observation model with fixed observation error
    }
  }
    "
      ,fill=TRUE)
  sink()
  closeAllConnections()
  
  sink("miss_ss.stan")
  
  cat("
 

/*----------------------- Functions --------------------------*/  
 
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int N; // Length of state and observation time series
    vector[N] y_miss; // Observations
    real z0; // Initial state value
    vector[N] light;
    vector[N] miss_vec; // index or location of missing values within the dataset
  }
  
/*----------------------- Parameters --------------------------*/
  /* Parameter block: defines the variables that will be sampled */
  parameters {
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real b0;
    real b1;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    vector[N] z; // State time series
  }
  

  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
  int count;
    // Prior distributions
    sdo ~ normal(0.1, 0.01);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    b1 ~ normal(0,5);
    
    // Distribution for the first state
   z[1] ~ normal(z0, sdp);
   //y_miss[1]~ normal(z[1],sdo);
   
    // Distributions for all other states
    count = 0;
    for(t in 2:N){
       z[t] ~ normal(b0+ z[t-1]*phi+light[t]*b1, sdp);
    if(miss_vec[t]==0){
      y_miss[t] ~ normal(z[t], sdo); // observation model with fixed observation error
    }
    }
  
  }
  
    "
      ,fill=TRUE)
  sink()
  closeAllConnections()
  
  
  
  MI<-function(m, y, light, ts){
    w<-as.data.frame(cbind(y,light,ts))
    data.amelia<-amelia(w, m = m, p2s=1, ts="ts", lags="y")
  }
  
  MI(5, y_full$x, light,y_full$time )
  
  
  #Assign NAs to the missing data in each list
  for(i in 1:m){
    z1<-data.amelia[[i]]$imputations$imp1$y_miss1
    z2<-data.amelia[[i]]$imputations$imp2$y_miss1
    z3<-data.amelia[[i]]$imputations$imp3$y_miss1
    z4<-data.amelia[[i]]$imputations$imp4$y_miss1
    z5<-data.amelia[[i]]$imputations$imp5$y_miss1
    q1<-z1[y_missing_integers[[i]]]
    q2<-z2[y_missing_integers[[i]]]
    q3<-z3[y_missing_integers[[i]]]
    q4<-z4[y_missing_integers[[i]]]
    q5<-z5[y_missing_integers[[i]]]
    q.ts<-ts[y_missing_integers[[i]]]
    q<-cbind(q1,q2,q3,q4,q5)
    q.mean<- apply(q[,-1], 1, mean)
    q.min<-apply(q[,-1], 1, min)
    q.max<-apply(q[,-1], 1, max)
    q.sd<- apply(q[,-1], 1, sd)
    x_obs_data<-as.data.frame(cbind(x,ts))
    x_imp_data<-as_tibble(cbind(q.mean,q.sd,q.ts, q.min,q.max))
    x_imp_data<-as.data.frame(x_imp_data %>% slice(1:length(y_missing_integers[[i]])))
    x_obs<-x[y_missing_integers[[i]]]
    x_combined<-as.data.frame(cbind(x_imp_data$q.mean,x_imp_data$q.sd,x_obs,x_imp_data$q.min,x_imp_data$q.max))
  }
  
  