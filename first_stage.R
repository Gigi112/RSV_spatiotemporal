# INITIALISE TO LOCAL CONTEXT
rm(list=ls())
# set work directory
setwd("~/OneDrive - Yale University/RSV/CAR/20000 iter/data")
# READ DATA
NYwide<- readRDS("NYwide.rds") # monthly time series of observed RSV hospitalizations from each ZIP code in NY
NJwide<- readRDS("NJwide.rds") # monthly time series of observed RSV hospitalizations from each ZIP code in NJ
CTwide<- readRDS("CTwide.rds") # monthly time series of observed RSV hospitalizations from each ZIP code in CT
nyunder5 <- readRDS("nyunder5.rds") #offset term, population under 5 of each ZIP code in NY
njunder5 <- readRDS("njunder5.rds") #offset term, population under 5 of each ZIP code in NJ
ctunder5 <- readRDS("ctunder5.rds") #offset term, population under 5 of each ZIP code in CT

#################################################################
#First stage model
#################################################################
firststage <- function(y,logpop,zip,date,n,m){
  library(rjags)
  library(coda)
  
model_string<-"
model{

for(i in 1:n.zip){
for(j in 1:n.date){
y[i,j]  ~ dpois(lambda[i,j])

log(lambda[i,j])<- logpop[i]+ beta0[i] + amp[i]*cos(2*pi*j/12 - phase[i]) + phi[i,j]

phi[i,j] ~ dnorm (0,tau3[i])
}

phase[i] <- (2*pi*(exp(trans[i]))) / (1+exp(trans[i]))
beta0[i] ~ dnorm(0, 0.001)
amp[i]<-exp(beta1[i])
beta1[i] ~ dnorm(0,0.001) 
trans[i] ~ dnorm(0,0.001)
tau3[i] ~ dgamma(1.00, 0.01)} 
}"

dataset <- list('y' = y,"logpop"=logpop, "n.zip"=zip,"n.date"=date,pi=pi)

model<-jags.model(textConnection(model_string),
                   data=dataset,
                   n.chains=2)
update(model, 
       n.iter=n)

modelresult<-coda.samples(model, variable.names=c("trans"),
                      thin=10,
                      n.iter=m)
result <- as.data.frame(modelresult[[1]])
first.result<<-result
}

#run model for CT
firststage(y=CTwide,logpop=ctunder5$log,zip=nrow(CTwide),date=ncol(CTwide),n=10000,m=50000)
cttransmean<- colMeans(first.result)  # posterior mean of phase estimates
var.ct <- c()
for (i in 1:275) {
  var.ct[i] <- var(first.result[,i])} # calculate variations of phase estimates
saveRDS(var.ct,"varct.rds")           # save variations for the second stage
saveRDS(cttransmean,"transct.rds")   # save phase estimates for the second stage

#run model for NY
firststage(y=NYwide,logpop=nyunder5$log,zip=nrow(NYwide),date=ncol(NYwide),n=10000,m=50000)
nytransmean<- colMeans(first.result)  # posterior mean of phase estimates
var.ny <- c()
for (i in 1:1745) {
  var.ny[i] <- var(first.result[,i])} # calculate variations of phase estimates
saveRDS(var.ny,"varny.rds")           # save variations for the second stage
saveRDS(first.result,"transny.rds")   # save phase estimates for the second stage

#run model for NJ
firststage(y=NJwide,logpop=njunder5$log,zip=nrow(NJwide),date=ncol(NJwide),n=10000,m=50000)
njtransmean<- colMeans(first.result)  # posterior mean of phase estimates
var.nj <- c()
for (i in 1:592) {
  var.nj[i] <- var(first.result[,i])} # calculate variations of phase estimates
saveRDS(var.nj,"varnj.rds")           # save variations for the second stage
saveRDS(first.result,"transnj.rds")   # save phase estimates for the second stage


