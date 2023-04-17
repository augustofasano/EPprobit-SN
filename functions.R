getParamsPFM = function(X , y ,nu2, moments = TRUE, predictive = FALSE, tolerance = 1e-2, maxIter = 1e4) {
  ######################################################
  # PRECOMPUTATION
  ######################################################
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  
  # compute H = X%*%V%*%t(X) and Omega_z directly or with Woodbury
  if(p<=n) {
    # define prior covariance matrix and its inverse
    Omega = diag(rep(nu2,p),p,p)
    invOmega = diag(rep(1/nu2,p),p,p)
    V = solve(crossprod(X)+invOmega)
    H = X%*%V%*%t(X)
    invOmZ = diag(1,nrow=n,ncol=n) - H # needed for ELBO
  } else{
    XXt = tcrossprod(X)
    invOmZ = solve(diag(1,nrow=n,ncol=n)+nu2*XXt) # needed for ELBO
    H = nu2*XXt%*%invOmZ
  }
  
  # compute optimal sigma2
  # h = diag(diag(H))
  sigma2 = as.double(1/(1-diag(H)), ncol = 1)
  sigma = sqrt(sigma2)
  
  # compute matrix to write the CAVI update in a vectorized form
  A = sigma2*H
  A[cbind(1:n,1:n)] = 0
  
  # other useful quantities needed for ELBO
  diagInvOmZ = diag(invOmZ)
  
  # initialization of variables
  mu = matrix(0,n,1)
  
  # inizialize coherently the vector of expectations meanZ
  musiRatio = as.double(mu/sigma)
  phiPhiRatio = exp(dnorm(musiRatio,log = T)-pnorm((2*y-1)*musiRatio,log.p = T))
  meanZ = mu + (2*y-1)*sigma*phiPhiRatio
  
  elbo = -Inf
  diff = 1
  nIter=0
  
  ######################################################
  # CAVI ALGORITHM
  ######################################################
  
  while(diff > tolerance & nIter < maxIter) {
    elboOld = elbo
    sumLogPhi = 0
    
    for(i in 1:n) {
      mu[i] = A[i,]%*%meanZ
      
      # compute first (needed for algorithm) and second (needed for ELBO) moments
      musiRatio = mu[i]/sigma[i]
      phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm((2*y[i]-1)*musiRatio, log.p = T))
      meanZ[i] = mu[i] + (2*y[i]-1)*sigma[i]*phiPhiRatio
      sumLogPhi = sumLogPhi + pnorm((2*y[i]-1)*musiRatio, log.p = T)
    }
    
    # computation of ELBO (up to an additive constant not depending on mu)
    elbo = -(t(meanZ)%*%invOmZ%*%meanZ -
               sum((meanZ^2)*diagInvOmZ))/2 -
      sum(meanZ*mu/sigma2) + sum((mu^2)/sigma2)/2 + sumLogPhi
    
    diff = abs(elbo-elboOld)
    nIter = nIter+1
    
    if(nIter%%100==0) {print(paste0("iter: ", nIter, ", ELBO: ", elbo))}
  }
  
  # get the optimal parameters of the normals before truncation, now that convergence has been reached
  mu = A%*%meanZ
  
  results = list(mu = mu, sigma2 = sigma2, nIter = nIter)
  
  ######################################################
  # (OPTIONAL) CLOSED-FORM MOMENTS' COMPUTATION
  ######################################################
  
  if(moments == TRUE) {
    # compute V and V%*%t(X), directly or with Woodbury
    if(p<=n) {
      diagV = diag(V) # V already computed
      VXt = V%*%t(X)
    } else{ # use Woodbury
      VXt = t(nu2*X)%*%invOmZ
      #V = Omega - VXt%*%(nu2*X)
      diagV = nu2*(1-colSums(t(VXt) * X))
    }
    
    musiRatio = mu/sigma
    phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm((2*y-1)*musiRatio, log.p = T))
    
    meanZ = mu + (2*y-1)*sigma*phiPhiRatio
    postVarZ = as.double(sigma2*(1-(2*y-1)*musiRatio*phiPhiRatio - phiPhiRatio^2))
    
    W = apply(VXt,1,function(x) sum(x*x*postVarZ))
    
    meanBeta = VXt%*%meanZ
    varBeta = diagV + W
    
    moments_PFM = list(meanBeta=meanBeta,varBeta=matrix(varBeta,ncol = 1))
    
    results = c(results,postMoments=moments_PFM)
  }
  
  if(predictive==T){
    # return also VXt and InvOmZ which are needed for the computation of predictive probabilities
    if(moments==F){
      # we weed to compute the quantities
      if(p<=n) {
        VXt = V%*%t(X)
      } else{
        VXt = t(nu2*X)%*%invOmZ
      }
    }
    # now we have quantities we need
    predQuant = list(VXt=VXt,invIXXt=invOmZ)
    results = c(results,predQuant=predQuant)
  }
  
  return(results)
}

rSUNpost = function(X,y,nu2,nSample) {
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  
  # get parameters useful for sampling
  Omega = diag(rep(nu2,p),p,p)
  invOmega = diag(rep(1/nu2,p),p,p)

  signX = X
  signX[y==0,] = -X[y==0,]
  Gamma_post_unnormalized = diag(1,n,n)+(nu2*signX)%*%t(signX)
  inv_Gamma_post_unnormalized = solve(Gamma_post_unnormalized)
  s = diag(sqrt(Gamma_post_unnormalized[cbind(1:n,1:n)]),n,n)
  s_1 = diag(1/s[cbind(1:n,1:n)],n,n)
  gamma_post = matrix(0,n,1) # because prior mean is set to 0
  Gamma_post = s_1%*%Gamma_post_unnormalized%*%s_1
  
  V = Omega-t(nu2*signX)%*%inv_Gamma_post_unnormalized%*%(signX*nu2)
  V = 0.5*(V+t(V))
  L = t(chol(V))
  
  # compute multiplicative coefficients for the truncated multivariate normal component
  coefTruncNorm = t(nu2*signX)%*%inv_Gamma_post_unnormalized%*%s
  
  # sample the multivariate normal component
  sampleMultNorm = matrix(rnorm(nSample*p),p,nSample)
  
  # sample the truncated multivariate normal component
  # if(n == 1) {
  #   sampleTruncNorm = matrix(rtruncnorm(n = nSample, a = -gamma_post, b = Inf, mean = 0, sd = 1), nrow = 1, ncol = nSample)
  # } else{
  #   sampleTruncNorm = mvrandn(l = -gamma_post, u = rep(Inf,n), Sig = Gamma_post, n = nSample)
  # } # this part is old now, after the update of TruncatedNormal to version 2.0
  sampleTruncNorm = t(rtmvnorm(n = nSample, mu = rep(0,n), sigma = Gamma_post, lb = -gamma_post, u = rep(Inf,n)))
  
  # combine the multivariate normal and truncated normal components
  sampleSUN = L%*%sampleMultNorm+coefTruncNorm%*%sampleTruncNorm
}

generateSyntheticDataNoCorr = function(nTrain, nTest, p, intercept = FALSE, beta = NULL, seed=1) {
  # generate n observations of p-variate regressors: if intercept == T the first column is all 1
  set.seed(seed)
  if(is.null(beta)){
    # randomly generate beta from uniform[-5,5]
    beta = matrix(runif(n=p, min = -5, max = 5), ncol = 1)
  }
  
  # generate train set
  if(intercept == T) {
    # generate train set and data and get rescaling quantities
    X_Train = matrix(rnorm((p-1)*nTrain,mean = 0, sd = 1), nrow = nTrain, ncol = p-1)
    meanTrain = apply(X_Train,2,mean)
    sdTrain = apply(X_Train,2,sd)
    
    # add intercept and standardise data
    X_Train = cbind(rep(1,nTrain),0.5*t(apply(X_Train,1, function(x) (x-meanTrain)/sdTrain)))
  } else{
    # generate train set and get rescaling quantities
    X_Train = matrix(rnorm(p*nTrain,mean = 0, sd = 1), nrow = nTrain, ncol = p)
    meanTrain = apply(X_Train,2,mean)
    sdTrain = apply(X_Train,2,sd)
    
    # standardise data
    X_Train = 0.5*t(apply(X_Train,1, function(x) (x-meanTrain)/sdTrain)) # so we have zero mean and std=0.5
  }
  
  # generate train data
  prob = pnorm(X_Train%*%beta)
  u = runif(nTrain,0,1)
  yTrain = as.integer(u <= prob)
  
  results = list(X_Train=X_Train,yTrain=yTrain,beta=beta)
  
  if(nTest>0){
    # generate test set
    if(intercept == T) {
      # generate test set
      X_Test = matrix(rnorm((p-1)*nTest,mean = 0, sd = 1), nrow = nTest, ncol = p-1)
      
      # add intercept and standardise data
      X_Test = cbind(rep(1,nTest),0.5*t(apply(X_Test,1, function(x) (x-meanTrain)/sdTrain)))
    } else{
      # generate test set
      X_Test = matrix(rnorm(p*nTest,mean = 0, sd = 1), nrow = nTest, ncol = p)
      
      # standardise data
      X_Test = 0.5*t(apply(X_Test,1, function(x) (x-meanTrain)/sdTrain)) # so we have zero mean and std=0.5
    }
    
    # generate test data
    prob = pnorm(X_Test%*%beta)
    u = runif(nTest,0,1)
    yTest = as.integer(u <= prob)
    
    results = append(results, list(X_Test=X_Test,yTest=yTest))
  }
  
  return(results)
}


# ------------------------------- EFFICIENT EP ---------------------------------
zeta1 = function(x){exp(dnorm(x, log = T) - pnorm(x,log.p = T))}
zeta2 = function(x,z1){-z1^2-x*z1}

getParamsEP = function(X,y,nu2,tolerance=1e-3,maxIter=1e3,fullVar=FALSE,predictive=FALSE,nPrint=100){
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  
  ### Initialization
  if(p<n){
    Q = diag(1/nu2,p,p)
    invQ = diag(nu2,p,p)
  }else{
    V = nu2*t(X)
  }
  
  r = rep(0,p)
  
  k = double(length = n)
  m = double(length = n)
  
  diff = 1
  nIter = 0
  
  ### Iterations
  
  if(p<n){
    
    while(diff > tolerance && nIter < maxIter){
      
      diff = 0.
      count = 0
      
      for(i in 1:n){
        
        xi = X[i,]
        
        Q_i = Q - k[i]*tcrossprod(xi)
        r_i = r - m[i]*xi
        
        Oxi = invQ%*%xi
        
        Oi = invQ + tcrossprod(Oxi)*k[i]/as.double(1.-k[i]*crossprod(xi,Oxi))
        
        Oixi = Oi%*%xi
        xiOixi = as.double(crossprod(xi,Oixi))
        
        if(xiOixi>0){
          
          r_iOixi = as.double(crossprod(r_i,Oixi))
          
          s = (2.*y[i]-1.)/sqrt(1.+xiOixi)
          tau = s*r_iOixi
          
          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)
          
          kNew = - z2/(1.+xiOixi+z2*xiOixi)
          mNew = s*z1 + kNew*r_iOixi + kNew*s*z1*xiOixi
          
          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff} #???
          
          k[i] = kNew
          m[i] = mNew
          
          r = r_i + m[i]*xi
          Q = Q_i + k[i]*tcrossprod(xi)
          
          invQ = Oi - tcrossprod(Oixi)*k[i]/(1.+k[i]*xiOixi)
        }else{
          count = count+1
          print(paste0(count," units skipped"))
        }
        
      }
      
      nIter = nIter + 1
      if(nIter %% nPrint == 0) {print(paste0("iteration ",nIter))}
    }
    
  } else {
    
    while(diff > tolerance && nIter < maxIter){
      
      diff = 0.
      count = 0
      
      for(i in 1:n){
        
        v = V[,i]
        xi = X[i,]
        xTv = crossprod(xi,v)[1]
        
        d = 1-k[i]*xTv
        w = v/d
        xTw = xTv/d
        
        if(xTw>0){
          
          r_iTw = crossprod(r,w) - m[i]*xTw
          
          s = (2*y[i]-1)/sqrt(1+xTw)
          tau = s*r_iTw
          
          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)
          
          kNew = as.double(-z2/(1 + xTw + z2*xTw))
          mNew = as.double(z1*s + kNew*r_iTw + kNew*z1*s*xTw)
          
          r = r + (mNew - m[i])*xi
          
          ratio = (k[i]-kNew)/(1.+(kNew-k[i])*xTv)
          V = V + (v*ratio)%*%crossprod(xi,V)
          
          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff}
          
          k[i] = kNew
          m[i] = mNew
          
        }else{
          count = count+1
          print(paste0(count," units skipped"))
        }
      }
      
      nIter = nIter + 1
      if(nIter %% nPrint == 0) {print(paste0("iteration ",nIter))}
    }
    
  }
  
  ### Posterior Approximate Moments
  
  if(p<n){
    
    meanBeta = invQ%*%r
    diagOmega = diag(invQ)
    
  }else{
    
    diagOmega = rep(nu2,p)*(rep(1,p) - rowSums(V*t(k*X)))
    meanBeta = nu2*(t(X)%*%m - V%*%(k*X%*%r))
  }
  
  results = list(meanBeta = meanBeta, diagOmega = diagOmega, 
                 nIter = nIter, kEP = k, mEP = m)
  
  if(fullVar==TRUE){
    if(p>=n) {
      invQ = nu2*(diag(1,p,p) - V%*%(k*X))
    }
    results = append(list(Omega=invQ),results)
  }
  
  if(predictive==TRUE){
    if(p>=n){
      results = append(list(V=V),results)
    } else{
      if(fullVar==FALSE){
        results = append(list(Omega=invQ),results)
      }
    }
  }
  
  return(results)
}

predictEP = function(paramsEP,xNew,nu2){
  p = length(paramsEP$meanBeta)
  n = length(paramsEP$kEP)
  
  if(p>=n){
    Xx = X%*%xNew
    KVtx = paramsEP$k*(t(paramsEP$V)%*%xNew)
    sd = as.double(sqrt(1+nu2*(sum(xNew^2)-sum(KVtx*Xx))))
  }else{
    sd = as.double(sqrt(1+t(xNew)%*%paramsEP$Omega%*%xNew))
  }
  predProb = as.double(pnorm(t(xNew)%*%paramsEP$meanBeta/sd))
  return(predProb)
}
