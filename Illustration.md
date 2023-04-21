# Introduction

As described in the [`README.md`](https://github.com/augustofasano/EPprobit-SN/blob/main/README.md) file, this tutorial contains general guidelines and code to **replicate the simulation studies** in the papers, computing approximate posterior moments and predictive probabilities via expectation propagation (EP).

# Preliminary operations

Once the file [`functions.R`](https://github.com/augustofasano/EPprobit-SN/blob/main/Functions.R) has been downloaded, set the working directory to the folder where it is located. Then, clean the workspace and load `it together with other useful packages.

```r
source("functions.R")
library("TruncatedNormal")
library("truncnorm")
library("EPGLM")
```

# Posterior inference for a fixed scenario n=100, p=800

We show how to replicate a specific scenario of the simulation studies.
Since the focus of the paper is on large p settings, **we consider the case n=100 and p=800**.
In the following, we also show how to replicate the whole simulation study, computing both approximate posterior moments and predictive probabilities in all the scenarios.

### Set useful quantities
We start our analysis by setting the hyperparameters of the model and defining the number of samples used for the i.i.d. sampler.

```r
nTrain = 100 # number of training observations that will be used to compute posterior moments
pnRatioValues = c(0.5,1,2,4,8)
pValues = nTrain*pnRatioValues
i=5 # consider the fifth scenario p=800
p = pValues[i]

nTest = 50 # number of test observations

nu2 = 25 # prior variance

nSample = 2e3 # number of samples for the i.i.d. sampler
```

We also allocate memory for the arrays where the predictive probabilities for the test units will be stored.

```r
predProbMC    = double(length = nTest)
predProbPFM   = double(length = nTest)
predProbEP    =  double(length = nTest)
predProbEPglm =  double(length = nTest)
```
Then we generate the data.

``` r
data = generateSyntheticDataNoCorr(nTrain = nTrain, nTest = nTest, p = p, intercept = T, seed = 100*i)
X = data$X_Train
y = data$yTrain
X_Test = data$X_Test
yTest = data$yTest
rm(data)
```

Now, we perform posterior inference according to the various methods considered in the papers and summarized in the [`README.md`](https://github.com/augustofasano/EPprobit-SN/blob/main/README.md) file.

### Posterior inference with the various methods

We start by getting 2000 samples from the **i.i.d. sampler** [(Durante, 2019)](https://academic.oup.com/biomet/article-abstract/106/4/765/5554418) via the function `rSUNpost`. From these samples, posterior means and standard deviations are computed via Monte Carlo approximation. The same is done also for posterior predictive probabilities.

````r
startTime = Sys.time()
betaSUN = rSUNpost(X=X,y=y,nu2=nu2,nSample=nSample)

# posterior moments
meanMC = apply(betaSUN,1,mean)
sdMC = apply(betaSUN,1,sd)
timeMC = difftime(Sys.time(), startTime, units=("secs"))[[1]]

# predictive probabilities
for(j in 1:nTest){
  xNew = matrix(X_Test[j,],ncol = 1)
  predProbMC[j] = mean(pnorm(t(xNew)%*%betaSUN))
}
timeMCpredProb = difftime(Sys.time(), startTime, units=("secs"))[[1]]
````

Secondly, we obtain the approximate posterior moments and predictive probabilities according to the **efficient expectation propagation (EP-EFF) implementation** presented in the papers. They are computed via the function `getParamsEP`.

````r
startTime = Sys.time()
tolerance = 1e-3 # tolerance to establish convergence
paramsEP = getParamsEP(X=X,y=y,nu2=nu2,tolerance=tolerance,predictive=T,fullVar = F,maxIter=1e4)
timeEP = difftime(Sys.time(), startTime, units=("secs"))[[1]]

# posterior moments
meanEP = paramsEP$meanBeta
sdEP = sqrt(paramsEP$diagOmega)

# predictive probabilities
for(j in 1:nTest){
  xNew = matrix(X_Test[j,],ncol = 1)
  predProbEP[j] = as.double(predictEP(paramsEP,xNew,nu2))
}
timeEPpredProb = difftime(Sys.time(), startTime, units=("secs"))[[1]]
````

Thirdly, we perform inference via the **partially factorized mean-field variational Bayes (PFM-VB)** approximation [(Fasano, Durante and Zanella, 2022)](https://academic.oup.com/biomet/article-abstract/109/4/901/6581071), implemented by the function `getParamsPFM`.

````r
tolerance = 1e-3 # tolerance to establish ELBO convergence
startTime = Sys.time()
paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,predictive=TRUE,tolerance=tolerance,maxIter=1e4)
timePFM = difftime(Sys.time(), startTime, units=("secs"))[[1]]

# posterior moments
meanPFM = paramsPFM$postMoments.meanBeta
sdPFM = sqrt(paramsPFM$postMoments.varBeta)

# predictive probabilities
# obtain a sample from q^*(z) to be used for the predictive probabilities of PFM
nSampleZ = 5e3
muTN = paramsPFM$mu
muTN[y==0] = -muTN[y==0]
sampleTruncNorm = matrix(rtruncnorm(nTrain*nSampleZ, a = 0, b = Inf, mean = muTN, sd = sqrt(paramsPFM$sigma2)), nrow = nTrain, ncol = nSampleZ, byrow = F )
sampleTruncNorm[y==0,] = -sampleTruncNorm[y==0,] # need to adjust the sign of the variables for which y_i is 0
  
for(j in 1:nTest){
  xNew = matrix(X_Test[j,],ncol = 1)
  Xx = X%*%xNew
  sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%paramsPFM$predQuant.invIXXt%*%Xx)))
    
  predProbPFM[j] = mean(pnorm((t(xNew)%*%paramsPFM$predQuant.VXt%*%sampleTruncNorm)/sd))
}
timePFMpredProb = difftime(Sys.time(), startTime, units=("secs"))[[1]]
````

Finally, we also compute the **EP approximate posterior moments and predictive probabilities according to routine implemented methods**, as the function `EPprobit` from the package `EPGLM`.
This clarifies to what extent the running time is decreased, without affecting the quality of the approximation, since the two implementations give the same approximate moments (up to numerical precision).

````r
startTime = Sys.time()
paramsEPglm = EPprobit(X=X,Y=y,s=nu2)
timeEPglm = difftime(Sys.time(), startTime, units=("secs"))[[1]]

# posterior moments
meanEPglm = paramsEPglm$m
sdEPglm = sqrt(diag(paramsEPglm$V))

# predictive probabilities
for(j in 1:nTest){
  xNew = matrix(X_Test[j,],ncol = 1)
  sd = as.double(sqrt(1+t(xNew)%*%paramsEPglm$V%*%xNew))
  predProbEPglm[j] = as.double(pnorm(t(xNew)%*%paramsEPglm$m/sd))
}
  
timeEPglm_predProb = difftime(Sys.time(), startTime, units=("secs"))[[1]]
````

### Comparison of posterior moments and running times
We can then compare the performances of the various methods.

We start by checking that our function `getParamsEP` obtains the same approximate moments as `EPprobit`.

```r
max(abs(meanEP-meanEPglm))
# [1] 3.264056e-14

max(abs(sdEP-sdEPglm))
# [1] 9.769963e-15

max(abs(predProbEP-predProbEPglm))
# [1] 1.665335e-15
```

Nevertheless, we see that the efficient implementation allows the computation of the approximate EP posterior moments at a massively lower computational time.

Time needed to compute posterior moments:

```r
round(timeEP,2)
round(timeEPglm,2)
round(timePFM,2)
```

Time needed to compute posterior moments and predictive probabilities:

```r
round(timeEPpredProb,2)
round(timeEPglm_predProb,2)
round(timePFMpredProb,2)
```


To check the quality of the EP and PFM-VB approximations, we compute the median absolute differences between the resulting approximate posterior means and standard deviations and the ones obtained with the i.i.d. sampler.

```r
round(quantile(abs(meanEP-meanMC),probs=0.5),2)
# 0.07
round(quantile(abs(meanPFM-meanMC),probs=0.5),2)
# 0.07

round(quantile(abs(sdEP-sdMC),probs=0.5),2)
# 0.05
round(quantile(abs(sdPFM-sdMC),probs=0.5),2)
# 0.05
```




# Code to reproduce the whole simulation study

We give here the code to reproduce the whole simulation study presented in the papers.

### Variables initialization

```r
samplesMC_acrossScens = list() # list of samples for the different scenarios
paramsPFM_acrossScens = list()
paramsEP_acrossScens = list()
paramsEPglm_acrossScens = list()

meanMC = list() # list of means for the parameters for different scenarios
meanPFM = list()
meanEP = list()
meanEPglm = list()

sdMC = list() # list of stds for the parameters for different scenarios
sdPFM = list()
sdEP = list()
sdEPglm = list()

predProbMC = matrix(nrow = nTest, ncol = length(pnRatioValues)) # matrix of predictive probabilities for 50 test obs for different scenarios
predProbPFM = matrix(nrow = nTest, ncol = length(pnRatioValues))
predProbEP =  matrix(nrow = nTest, ncol = length(pnRatioValues))
predProbEPglm =  matrix(nrow = nTest, ncol = length(pnRatioValues))

timeMC= double(length = length(pValues))
timePFM = double(length = length(pValues))
timeEP = double(length = length(pValues))
timeEPglm = double(length = length(pValues))

timeMCpredProb = double(length = length(pValues))
timePFMpredProb = double(length = length(pValues))
timeEPpredProb = double(length = length(pValues))
timeEPglm_predProb = double(length = length(pValues))
```

### Run the simulations and save the output

```r
for(i in 1:length(pValues)) {
  ######################################################
  # GENERATE  THE DATA
  ######################################################
  p = pValues[i]
  print(paste(" p:",p))
  data = generateSyntheticDataNoCorr(nTrain = nTrain, nTest = nTest, p = p, intercept = T, seed = 100*i)
  X = data$X_Train
  y = data$yTrain
  X_Test = data$X_Test
  yTest = data$yTest
  rm(data)
  
  ######################################################
  # GET MC SAMPLES
  ######################################################
  # get posterior samples
  startTime = Sys.time()
  betaSUN = rSUNpost(X=X,y=y,nu2=nu2,nSample=nSample)
  meanMC[[i]] = apply(betaSUN,1,mean)
  sdMC[[i]] = apply(betaSUN,1,sd)
  timeMC[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # predictive probabilities
  for(j in 1:nTest){
    xNew = matrix(X_Test[j,],ncol = 1)
    
    predProbMC[j,i] = mean(pnorm(t(xNew)%*%betaSUN))
  }
  timeMCpredProb[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # store quantities
  samplesMC_acrossScens[[i]] = betaSUN
  
  ######################################################
  # GET PFM POSTERIOR MOMENTS
  ######################################################
  tolerance = 1e-3 # tolerance to establish ELBO convergence
  # get optimal parameters and moments and predictive probabilities
  startTime = Sys.time()
  paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,predictive=TRUE,tolerance=tolerance,maxIter=1e4)
  timePFM[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # predictive probabilities
  # obtain a sample from q^*(z) to be used for the predictive probabilities of PFM
  nSampleZ = 2e3# 5e3
  muTN = paramsPFM$mu
  muTN[y==0] = -muTN[y==0]
  sampleTruncNorm = matrix(rtruncnorm(nTrain*nSampleZ, a = 0, b = Inf, mean = muTN, sd = sqrt(paramsPFM$sigma2)), nrow = nTrain, ncol = nSampleZ, byrow = F )
  sampleTruncNorm[y==0,] = -sampleTruncNorm[y==0,] # need to adjust the sign of the variables for which y_i is 0
  
  for(j in 1:nTest){
    xNew = matrix(X_Test[j,],ncol = 1)
    Xx = X%*%xNew
    sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%paramsPFM$predQuant.invIXXt%*%Xx)))
    
    predProbPFM[j,i] = mean(pnorm((t(xNew)%*%paramsPFM$predQuant.VXt%*%sampleTruncNorm)/sd))
  }
  timePFMpredProb[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # store quantities
  meanPFM[[i]] = paramsPFM$postMoments.meanBeta
  sdPFM[[i]] = sqrt(paramsPFM$postMoments.varBeta)
  
  ######################################################
  # GET EP POSTERIOR MOMENTS
  ######################################################
  startTime = Sys.time()
  tolerance = 1e-3 # tolerance to establish ELBO convergence
  paramsEP = getParamsEP(X=X,y=y,nu2=nu2,tolerance=tolerance,predictive=T,fullVar = F,maxIter=1e4)
  timeEP[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # predictive probabilities
  for(j in 1:nTest){
    xNew = matrix(X_Test[j,],ncol = 1)
    predProbEP[j,i] = as.double(predictEP(paramsEP,xNew,nu2))
  }
  
  timeEPpredProb[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # store quantities
  meanEP[[i]] = paramsEP$meanBeta
  sdEP[[i]] = sqrt(paramsEP$diagOmega)
  
  
  ######################################################
  # GET EP POSTERIOR MOMENTS EPGLM
  ######################################################
  startTime = Sys.time()
  tolerance = 1e-3 # tolerance to establish ELBO convergence
  paramsEPglm = EPprobit(X=X,Y=y,s=nu2)
  timeEPglm[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # predictive probabilities
  for(j in 1:nTest){
    xNew = matrix(X_Test[j,],ncol = 1)
    sd = as.double(sqrt(1+t(xNew)%*%paramsEPglm$V%*%xNew))
    predProbEPglm[j,i] = as.double(pnorm(t(xNew)%*%paramsEPglm$m/sd))
  }
  
  timeEPglm_predProb[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  
  # store quantities
  meanEPglm[[i]] = paramsEPglm$m
  sdEPglm[[i]] = sqrt(diag(paramsEPglm$V))
}

save(meanMC,meanPFM,meanEP,meanEPglm,
     predProbMC,predProbPFM,predProbEP,predProbEPglm,
     sdMC,sdPFM,sdEP,sdEPglm,
     samplesMC_acrossScens,
     timeMC,timePFM,timeEP,timeEPglm,
     timeMCpredProb,timePFMpredProb,timeEPpredProb,timeEPglm_predProb,
     file = "outputSimStudy.RData")
```
### Reproduce the results of Fasano et al., Book of short papers - SIS 2023

#### Running times to compute posterior moments (Table 1)

```r
rm(list=ls())
load("outputSimStudy.RData")

round(timeEP,3)
# [1] 0.110 0.016 0.031 0.050 0.087
round(timeEPglm,3)
# [1]   0.067   0.416   3.185  24.363 140.244
round(timePFM,3)
# [1] 0.110 0.063 0.005 0.007 0.011
```

#### Get absolute differences between i.i.d. posterior moments and the ones obtained via EP and PFM-VB (Figure 1)

First, we need to create an appropriate data frame.

```r
rm(list=ls())

load("outputSimStudy.RData")
nScen = length(meanMC)

meanAbsDiffPFM=matrix(0,nScen,3)
meanAbsDiffEP=matrix(0,nScen,3)
sdAbsDiffPFM=matrix(0,nScen,3)
sdAbsDiffEP=matrix(0,nScen,3)

for (i in 1:nScen) {  
  meanAbsDiffPFM[i,] = quantile(abs(meanPFM[[i]]-meanMC[[i]]),probs=c(0.25,0.5,0.75))
  meanAbsDiffEP[i,] = quantile(abs(meanEP[[i]]-meanMC[[i]]),probs=c(0.25,0.5,0.75))
  
  sdAbsDiffPFM[i,] = quantile(abs(sdPFM[[i]]-sdMC[[i]]),probs=c(0.25,0.5,0.75))
  sdAbsDiffEP[i,] = quantile(abs(sdEP[[i]]-sdMC[[i]]),probs=c(0.25,0.5,0.75))
}

Plot_dataset<-data.frame(c(rep("Mean",5),rep("Standard deviation",5),rep("Mean",5),rep("Standard deviation",5)))
Plot_dataset<-cbind(Plot_dataset,c(rep("PFM-VB",10),rep("EP-EFF",10)),rep(c(50,100,200,400,800),4))

Plot_dataset<-cbind(Plot_dataset,rep(0,10*2))
Plot_dataset<-cbind(Plot_dataset,rep(0,10*2))
Plot_dataset<-cbind(Plot_dataset,rep(0,10*2))

for (i in 1:5){
  Plot_dataset[i,4]<-	meanAbsDiffPFM[i,1]
  Plot_dataset[i,5]<-	meanAbsDiffPFM[i,2]
  Plot_dataset[i,6]<-	meanAbsDiffPFM[i,3]
}

for (i in 6:10){
  Plot_dataset[i,4]<-	sdAbsDiffPFM[i-5,1]
  Plot_dataset[i,5]<-	sdAbsDiffPFM[i-5,2]
  Plot_dataset[i,6]<-	sdAbsDiffPFM[i-5,3]
}

for (i in 11:15){
  Plot_dataset[i,4]<-	meanAbsDiffEP[i-10,1]
  Plot_dataset[i,5]<-	meanAbsDiffEP[i-10,2]
  Plot_dataset[i,6]<-	meanAbsDiffEP[i-10,3]
}

for (i in 16:20){
  Plot_dataset[i,4]<-	sdAbsDiffEP[i-15,1]
  Plot_dataset[i,5]<-	sdAbsDiffEP[i-15,2]
  Plot_dataset[i,6]<-	sdAbsDiffEP[i-15,3]
}

colnames(Plot_dataset)<-c("moment","method","p","low","median","high")

Plot_dataset$method<-factor(Plot_dataset$method,levels=c("EP-EFF","PFM-VB"))

Plot_dataset$moment<-factor(Plot_dataset$moment,levels=c("Mean","Standard deviation"))
```
And then we plot the quantities.

```r
library(ggplot2)
library(reshape)
library(scales)

ggplot(data=Plot_dataset)+
  geom_ribbon(aes(x=p,ymin=low, ymax=high, group=method),alpha=0.05)+
  geom_line(aes(x=p, y=median,lty=method,color=method),lwd=0.5)+
  scale_color_manual(values=rep("#4a4a4a",2))+
  scale_linetype_manual(values=rep("solid",2))+
  geom_point(aes(x=p, y=median,shape=method,color=method),size=2)+
  scale_shape_manual(values=c(16,17))+
  facet_wrap(moment~.)+
  theme_bw()+xlab("p")+ylab("")+
  theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),axis.text.y =element_text(size=7),plot.margin = margin(0.1, 0.1, 0.05, -0.2, "cm"),legend.position = "bottom",legend.title=element_blank())
```

![alt text](https://github.com/augustofasano/EPprobit-SN/blob/main/fig1.png)


### Reproduce the results of Fasano et al., Book of short papers - CLADAG 2023

#### Running times to run the algorithm and compute the approximate posterior predictive probabilities

```r
rm(list=ls())
load("outputSimStudy.RData")

round(timeEPpredProb,2)
# [1] 0.12 0.02 0.03 0.06 0.10
round(timeEPglm_predProb,2)
# [1]   0.07   0.42   3.19  24.37 140.29
round(timePFMpredProb,2)
# [1] 0.23 0.19 0.13 0.14 0.14
```

#### Get absolute differences between posterior predictive probabilities arising from i.i.d. posterior samples and the ones obtained via EP and PFM-VB (Figure 1)

First, we need to create an appropriate data frame.

```r
rm(list=ls())

load("outputSimStudy.RData")
nScen = length(meanMC)

predProbAbsDiffPFM=matrix(0,nScen,3)
predProbAbsDiffEP=matrix(0,nScen,3)

for (i in 1:nScen) {
  predProbAbsDiffPFM[i,] = quantile(abs(predProbPFM[,i]-predProbMC[,i]),probs=c(0.25,0.5,0.75))
  predProbAbsDiffEP[i,] = quantile(abs(predProbEP[,i]-predProbMC[,i]),probs=c(0.25,0.5,0.75))
}

PlotPred_dataset <-data.frame(c(rep("Predictive probabilities",10)))
PlotPred_dataset <-cbind(PlotPred_dataset,c(rep("PFM-VB",5),rep("EP-EFF",5)),rep(c(50,100,200,400,800),2))
for (i in 1:5){
  PlotPred_dataset[i,4]<-	predProbAbsDiffPFM[i,1]
  PlotPred_dataset[i,5]<-	predProbAbsDiffPFM[i,2]
  PlotPred_dataset[i,6]<-	predProbAbsDiffPFM[i,3]
}

for (i in 6:10){
  PlotPred_dataset[i,4]<-	predProbAbsDiffEP[i-5,1]
  PlotPred_dataset[i,5]<-	predProbAbsDiffEP[i-5,2]
  PlotPred_dataset[i,6]<-	predProbAbsDiffEP[i-5,3]
}


colnames(PlotPred_dataset)<-c("quantity","method","p","low","median","high")

PlotPred_dataset$method<-factor(PlotPred_dataset$method,levels=c("EP-EFF","PFM-VB"))

PlotPred_dataset$quantity<-factor(PlotPred_dataset$quantity,levels=c("Predictive probabilities"))
```
And then we plot the quantities.

```r
ggplot(data=PlotPred_dataset)+
  geom_ribbon(aes(x=p,ymin=low, ymax=high, group=method),alpha=0.05)+
  geom_line(aes(x=p, y=median,lty=method,color=method),lwd=0.5)+
  scale_color_manual(values=rep("#4a4a4a",2))+
  scale_linetype_manual(values=rep("solid",2))+
  geom_point(aes(x=p, y=median,shape=method,color=method),size=2)+
  scale_shape_manual(values=c(16,17))+
  facet_wrap(quantity~.)+
  theme_bw()+xlab("p")+ylab("")+
  theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),axis.text.y =element_text(size=7),plot.margin = margin(0.1, 0.1, 0.05, -0.2, "cm"),legend.position = "bottom",legend.title=element_blank())
```

![alt text](https://github.com/augustofasano/EPprobit-SN/blob/main/figPredProbs.png)



