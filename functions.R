getParamsPFM = function(X,y,nu2,moments = TRUE,tolerance = 1e-2, maxIter = 1e4) {
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
  n = nTrain + nTest
  if(is.null(beta)){
    # randomly generate beta from uniform[-5,5]
    beta = matrix(runif(n=p, min = -5, max = 5), ncol = 1)
  }
  
  if(intercept == T) {
    X = matrix(rnorm((p-1)*n,mean = 0, sd = 1), nrow = n, ncol = p-1)
    meanTrain = apply(X[1:nTrain,],2,mean)
    sdTrain = apply(X[1:nTrain,],2,sd)
    X = cbind(rep(1,n),0.5*t(apply(X,1, function(x) (x-meanTrain)/sdTrain)))
  } else{
    X = matrix(rnorm(p*n,mean = 0, sd = 1), nrow = n, ncol = p)
    meanTrain = apply(X[1:nTrain,],2,mean)
    sdTrain = apply(X[1:nTrain,],2,sd)
    X = 0.5*t(apply(X,1, function(x) (x-meanTrain)/sdTrain)) # so we have zero mean and std=0.5
  }
  
  prob = pnorm(X%*%beta)
  u = runif(n,0,1)
  y = as.integer(u <= prob)
  
  return(list(X=X,y=y,beta=beta))
}

generateSyntheticDataNoCorrNoStd = function(nTrain, nTest, p, intercept = FALSE, beta = NULL, seed=1) {
  # generate n observations of p-variate regressors: if intercept == T the first column is all 1
  set.seed(seed)
  n = nTrain + nTest
  if(is.null(beta)){
    # randomly generate beta from uniform[-5,5]
    beta = matrix(runif(n=p, min = -5, max = 5), ncol = 1)
  }
  
  if(intercept == T) {
    X = matrix(rnorm((p-1)*n,mean = 0, sd = 1), nrow = n, ncol = p-1)
    # meanTrain = apply(X[1:nTrain,],2,mean)
    # sdTrain = apply(X[1:nTrain,],2,sd)
    # X = cbind(rep(1,n),0.5*t(apply(X,1, function(x) (x-meanTrain)/sdTrain)))
    X = cbind(rep(1,n),X)
  } else{
    X = matrix(rnorm(p*n,mean = 0, sd = 1), nrow = n, ncol = p)
    # meanTrain = apply(X[1:nTrain,],2,mean)
    # sdTrain = apply(X[1:nTrain,],2,sd)
    # X = 0.5*t(apply(X,1, function(x) (x-meanTrain)/sdTrain)) # so we have zero mean and std=0.5
  }
  
  prob = pnorm(X%*%beta)
  u = runif(n,0,1)
  y = as.integer(u <= prob)
  
  return(list(X=X,y=y,beta=beta))
}


# ------------------------------- EFFICIENT EP ---------------------------------
zeta1 = function(x){exp(dnorm(x, log = T) - pnorm(x,log.p = T))}
zeta2 = function(x,z1){-z1^2-x*z1}

getParamsEP = function(X,y,nu2,tolerance=1e-3,maxIter=1e3,fullVar=FALSE,nPrint=100){
  
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
      Omega = nu2*(diag(1,p,p) - V%*%(k*X))
    }
    results = append(list(Omega=Omega),results)
  }
  
  return(results)
}

getParamsEP_bak = function(X1,y1,X0,y0,om2,zT,tolerance=1e-3,maxIter=1e3,
                       fullVar=FALSE,predictive=FALSE,nPrint=1000){
  
  n1 = dim(X1)[1]
  n0 = dim(X0)[1]
  p = dim(X0)[2]
  
  ### Pre-Computations
  
  r = crossprod(X1,y1) + c(-zT/om2,rep(0,p-1))
  if(p<n1){
    Omega0 = solve(diag(1./om2,p,p) + crossprod(X1))
    beta0 = Omega0%*%r
  }else{
    Lambda = solve(diag(1.,n1,n1)+om2*tcrossprod(X1))
    Omega0 = diag(om2,p,p) - om2*crossprod(X1,Lambda%*%X1)*om2
    beta0 = om2*(r - crossprod(X1,Lambda%*%(X1%*%r))*om2)
  }
  # Omega0 = 0.5*(Omega0+t(Omega0))
  
  if(p<n0){
    Omega = Omega0
    if(p<n1){
      Q = solve(Omega)
    }else{
      Q = diag(1./om2,p,p) + crossprod(X1)
    }
    logDetQ0 = determinant(Q, logarithm = TRUE)
  }else{
    if(p<n1){
      U = tcrossprod(Omega0,X0)
    } else {
      U = om2*t(X0) - om2*crossprod(X1,Lambda%*%tcrossprod(X1,X0))*om2
    }
    U0 = U
  }
  
  logZ0 = 0.5*crossprod(r,beta0)
  
  ### Initialization
  
  logZ = double(length = n0)
  k = double(length = n0)
  m = double(length = n0)
  
  diff = 1
  nIter = 0
  
  ### Iterations
  
  if(p<n0){
    
    while(diff > tolerance && nIter < maxIter){
      
      diff = 0.
      count = 0
      logZ = 0.
      
      for(i in c(1:n0)){
        
        xi = X0[i,]
        
        r_i = r - m[i]*xi
        Q_i = Q - k[i]*tcrossprod(xi)
        
        Oxi = Omega%*%xi
        
        Oi = Omega + tcrossprod(Oxi)*k[i]/as.double(1.-k[i]*crossprod(xi,Oxi))
        
        Oixi = Oi%*%xi
        xiOixi = as.double(crossprod(xi,Oixi))
        
        if(xiOixi>0){
          
          r_iOixi = as.double(crossprod(r_i,Oixi))
          
          s = (2.*y0[i]-1.)/sqrt(1.+xiOixi)
          tau = s*r_iOixi
          
          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)
          
          kNew = - z2/(1.+xiOixi+z2*xiOixi)
          mNew = s*z1 + kNew*r_iOixi + kNew*s*z1*xiOixi
          
          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff}
          
          k[i] = kNew
          m[i] = mNew
          
          logZ[i] = pnorm(tau,log.p = T) + 0.5*log(1.+k[i]*xiOixi) +
            0.5*((r_iOixi)^2)/xiOixi - 0.5*((m[i] + r_iOixi/xiOixi)^2)*xiOixi/(1.+k[i]*xiOixi)
          
          r = r_i + m[i]*xi
          Q = Q_i + k[i]*tcrossprod(xi)
          
          Omega = Oi - tcrossprod(Oixi)*k[i]/(1.+k[i]*xiOixi)
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
      logZ = 0.
      
      for(i in c(1:n0)){
        
        u = U[,i]
        xi = X0[i,]
        xTu = crossprod(xi,u)[1]
        
        d = 1-k[i]*xTu
        w = u/d
        xTw = xTu/d
        
        if(xTw>0){
          
          r_iTw = crossprod(r,w) - m[i]*xTw
          
          s = (2*y0[i]-1)/sqrt(1+xTw)
          tau = s*r_iTw
          
          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)
          
          kNew = as.double(-z2/(1 + xTw + z2*xTw))
          mNew = as.double(z1*s + kNew*r_iTw + kNew*z1*s*xTw)
          
          r = r + (mNew - m[i])*xi
          
          ratio = (k[i]-kNew)/(1.+(kNew-k[i])*xTu)
          U = U + (u*ratio)%*%crossprod(xi,U)
          
          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff}
          
          k[i] = kNew
          m[i] = mNew
          
          logZ[i] = pnorm(tau,log.p = T) + 0.5*log(1.+k[i]*xTw) +
            0.5*((r_iTw)^2)/xTw - 0.5*((m[i]+r_iTw/xTw)^2)*xTw/(1.+k[i]*xTw)
          
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
  
  if(p<n0){
    
    meanBeta = Omega%*%r
    diagOmega = diag(Omega)
    
    logDetQ = determinant(Q, logarithm = TRUE)
    logML = sum(logZ) + 0.5*crossprod(r,meanBeta) - logZ0 + 
      0.5*logDetQ0$modulus[1] - 0.5*logDetQ$modulus[1] 
    
  }else{
    
    diagOmega = diag(Omega0) - rowSums(U*t(k*t(U0)))
    meanBeta = beta0 + U0%*%m - U%*%(k*crossprod(U0,r))
    
    logDet = determinant(diag(1,n0,n0) + k*X0%*%U0, logarithm = TRUE)
    logML = sum(logZ) + 0.5*crossprod(r,meanBeta) - logZ0 - 0.5*logDet$modulus[1]
  }
  
  results = list(meanBeta = meanBeta, diagOmega = diagOmega, logML = logML, 
                 nIter = nIter, kEP = k, mEP = m)
  
  if(fullVar==TRUE){
    if(p>=n0) {
      Omega = Omega0 - U%*%(k*t(U0))
    }
    results = append(list(Omega=Omega),results)
  }
  
  if(predictive==TRUE){
    if(p>=n0) {
      results = c(results,postPredictive=list(U=U,U0=U0))
    }
  }
  
  return(results)
}
