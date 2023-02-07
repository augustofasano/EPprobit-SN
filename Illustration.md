# Introduction

As described in the [`README.md`](https://github.com/augustofasano/EPprobit-SN/blob/main/README.md) file, this tutorial contains general guidelines and code to **replicate the simulation study** in the paper.

# Preliminary operations

Once the file [`functions.R`](https://github.com/augustofasano/EPprobit-SN/blob/main/Functions.R) has been downloaded, set the working directory to the folder where it is located. Then, clean the workspace and load `it together with other useful packages.

```r
source("functions.R")
library("TruncatedNormal")
library("EPGLM")
```

# Posterior inference for a fixed scenario n=100, p=800

We show how to replicate a specific scenario of the simulation study.
Since the focus of the paper is on large p settings, we consider the case n=100 and p=800.
In the following we also show how to replicate the whole simulation study with all the scenarios.

### Set useful quantities
We start our analysis setting the hyperparameters of the model and defining the number of samples used for the i.i.d. sampler.

```r
n = 100
pnRatioValues = c(0.5,1,2,4,8)
pValues = n*pnRatioValues
i=5 # consider the fifth scenario p=800
p = pValues[i]

nu2 = 25 # prior variance

nSample = 2e3 # number of samples for the i.i.d. sampler
```
Then we generate the data.

``` r
data = generateSyntheticDataNoCorr(nTrain = n, nTest = 0, p = p, intercept = T, seed = 100*i)
X = data$X
y = data$y
rm(data)
```

Now, we perform posterior inference according to the various methods considered in the paper and summarized in the [`README.md`](https://github.com/augustofasano/EPprobit-SN/blob/main/README.md) file.

### Posterior inference with the various methods

We start by getting 2000 samples from the **i.i.d. sampler** [(Durante, 2019)](https://academic.oup.com/biomet/article-abstract/106/4/765/5554418) via the function `rSUNpost`. From these samples, posterior means and standard deviations are computed via Monte Carlo approximation.

````r
startTime = Sys.time()
betaSUN = rSUNpost(X=X,y=y,nu2=nu2,nSample=nSample)
timeMC = difftime(Sys.time(), startTime, units=("secs"))[[1]]

meanMC = apply(betaSUN,1,mean)
sdMC = apply(betaSUN,1,sd)
````

Secondly, we obtain the approximate posterior moments according the the **efficient expectation propagation (EP) implementation** presented in Algorithms 1 and 2 in the paper. They are computed via the function `getParamsEP`.

````r
startTime = Sys.time()
tolerance = 1e-3 # tolerance to establish convergence
paramsEP = getParamsEP(X=X,y=y,nu2=nu2,tolerance=tolerance,maxIter=1e4)
timeEP = difftime(Sys.time(), startTime, units=("secs"))[[1]]

meanEP = paramsEP$meanBeta
sdEP = sqrt(paramsEP$diagOmega)
````


Thirdly, we perform inference via the **partially factorized mean-field variational Bayes (PFM-VB)** approximation [(Fasano, Durante and Zanella, 2022)](https://academic.oup.com/biomet/article-abstract/109/4/901/6581071), implemented by the function `getParamsPFM`.

````r
tolerance = 1e-3 # tolerance to establish ELBO convergence
startTime = Sys.time()
paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,tolerance=tolerance,maxIter=1e4)
timePFM = difftime(Sys.time(), startTime, units=("secs"))[[1]]

meanPFM = paramsPFM$postMoments.meanBeta
sdPFM = sqrt(paramsPFM$postMoments.varBeta)
````

Finally, we also compute the **EP approximate moments according to routine implemented methods**, as the function `EPprobit` from the package `EPGLM`.
This clarifies to what extent the running time is decreased, without affecting the quality of the approximation, since the two implementations give the same approximate moments (up to numerical precision).

````r
startTime = Sys.time()
paramsEPglm = EPprobit(X=X,Y=y,s=nu2)
timeEPglm = difftime(Sys.time(), startTime, units=("secs"))[[1]]

meanEPglm = paramsEPglm$m
sdEPglm = sqrt(diag(paramsEPglm$V))
````

### Comparison of posterior moments and running times
We can then compare the performances of the various methods.

We start by checking that our function `getParamsEP` obtains the same approximate moments as `EPprobit`.

```r
max(abs(meanEP-meanEPglm))
# [1] 3.264056e-14

max(abs(sdEP-sdEPglm))
# [1] 9.769963e-15
```

Nevertheless, we see that the efficient implementation allows the computation of the approximate EP posterior moments at massively lower computational time.


```r
round(timeEP,2)
# [1] 0.09
round(timeEPglm,2)
# [1] 161.14
round(timePFM,2)
# [1] 0.01
```

To check the quality of the EP and PFM-VB approximations, we compute the median absolute differences between the resulting approximate posterior means and standard deviations and the ones obtained with the i.i.d. sampler.

```r
quantile(abs(meanEP-meanMC),probs=0.5)
# 0.07177436 
quantile(abs(meanPFM-meanMC),probs=0.5)
# 0.07224323

quantile(abs(sdEP-sdMC),probs=0.5)
# 0.05017391
quantile(abs(sdPFM-sdMC),probs=0.5)
# 0.04931011
```




# Code to reproduce the whole simulation study

We give here the code to reproduce the whole simulation study in the paper.

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


timeMC= double(length = length(pValues))
timePFM = double(length = length(pValues))
timeEP = double(length = length(pValues))
timeEPglm = double(length = length(pValues))
```

### Run the simulations and save the output

```r
for(i in 1:length(pValues)) {
  ######################################################
  # GENERATE  THE DATA
  ######################################################
  p = pValues[i]
  print(paste(" p:",p))
  data = generateSyntheticDataNoCorr(nTrain = n, nTest = 0, p = p, intercept = T, seed = 100*i)
  X = data$X
  y = data$y
  rm(data)
  
  ######################################################
  # GET MC SAMPLES
  ######################################################
  startTime = Sys.time()
  betaSUN = rSUNpost(X=X,y=y,nu2=nu2,nSample=nSample)
  timeMC[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  samplesMC_acrossScens[[i]] = betaSUN
  meanMC[[i]] = apply(betaSUN,1,mean)
  sdMC[[i]] = apply(betaSUN,1,sd)
  
  ######################################################
  # GET PFM POSTERIOR MOMENTS
  ######################################################
  tolerance = 1e-3 # tolerance to establish ELBO convergence
  # get optimal parameters and moments and predictive probabilities
  startTime = Sys.time()
  paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,tolerance=tolerance,maxIter=1e4)
  timePFM[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  meanPFM[[i]] = paramsPFM$postMoments.meanBeta
  sdPFM[[i]] = sqrt(paramsPFM$postMoments.varBeta)
  
  ######################################################
  # GET EP POSTERIOR MOMENTS
  ######################################################
  startTime = Sys.time()
  tolerance = 1e-3 # tolerance to establish ELBO convergence
  paramsEP = getParamsEP(X=X,y=y,nu2=nu2,tolerance=tolerance,maxIter=1e4)
  timeEP[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  meanEP[[i]] = paramsEP$meanBeta
  sdEP[[i]] = sqrt(paramsEP$diagOmega)
  
  
  ######################################################
  # GET EP POSTERIOR MOMENTS EPGLM
  ######################################################
  startTime = Sys.time()
  tolerance = 1e-3 # tolerance to establish ELBO convergence
  paramsEPglm = EPprobit(X=X,Y=y,s=nu2)
  timeEPglm[i] = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  meanEPglm[[i]] = paramsEPglm$m
  sdEPglm[[i]] = sqrt(diag(paramsEPglm$V))
}

save(meanMC,meanPFM,meanEP,meanEPglm,
     sdMC,sdPFM,sdEP,sdEPglm,
     samplesMC_acrossScens,
     timeMC,timePFM,timeEP,timeEPglm,
     file = "outputSimStudy.RData")
```

### Running times (Table 1)

```r
rm(list=ls())
load("outputSimStudy.RData")

round(timeEP,2)
# [1] 0.11 0.02 0.03 0.05 0.09
round(timeEPglm,2)
# [1]   0.07   0.43   3.09  25.32 161.14
round(timePFM,2)
# [1] 0.10 0.06 0.01 0.01 0.01
```

### Get absolute differences between i.i.d. posterior moments and the ones obtained via EP and PFM-VB (Figure 1)
First we need to create an appropriate data frame.

```r
Plot_dataset<-data.frame(c(rep("Mean",5),rep("Standard deviation",5),rep("Mean",5),rep("Standard deviation",5)))
Plot_dataset<-cbind(Plot_dataset,c(rep("PFM-VB",10),rep("EP",10)),rep(c(50,100,200,400,800),4))

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

Plot_dataset$method<-factor(Plot_dataset$method,levels=c("EP","PFM-VB"))

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
  scale_color_manual(values=c(rep("#4a4a4a",2),rep("#9e9d9d",2)))+
  scale_linetype_manual(values=c("solid","solid",rep("longdash",2)))+
  geom_point(aes(x=p, y=median,shape=method,color=method),size=2)+
  scale_shape_manual(values=c(16,17))+
  facet_wrap(moment~.)+
  theme_bw()+xlab("p")+ylab("")+
  theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),axis.text.y =element_text(size=7),plot.margin = margin(0.1, 0.1, 0.05, -0.2, "cm"),legend.position = "bottom",legend.title=element_blank())
```

![alt text](https://github.com/augustofasano//EPprobit-SN/blob/main/fig1.png)
