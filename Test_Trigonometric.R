rm(list=ls()) 



####### 




#### Loading packages

library(Matrix)
library(nloptr)
library(minqa)
library(optimx)
library(rootSolve)
library(nlme)
library(mvtnorm)
library(MASS)

final_res = list()

for (cond in 1:3) {
  


#############

res_mat = matrix(0,50,8)
for (run in 1:2) {
  

set.seed(run)
n=3;  #number of profiles

fun=function(x,i) { #defining all the profiles except target
  if     (i==1){
    
    return(25+2*cos(x) +rnorm(length(x),0,1)) }
  else if(i==2){
    
    return( -25+2* sin(x) +rnorm(length(x),0,1))    }
  
}
fun2=function(x){ #target profile
  return(10+2 * exp(-0.2*x) * cos(x)+rnorm(length(x),0,1))
}

num=c(5,4,6) #vector representing number of repetition for each profile
rep_factor =c() # this vector is used to transform off diagonal element denominators from sqrt(m_i*m_n) to max(m_i,m_n)(refer to equation 13)
for (i in 1:n-1) {
  rep_factor[i] = sqrt(num[i]*num[n])/max(num[i] , num[n])
}
trains <- seq(0,6*pi,length.out=45) #time points
trainy1=lapply(1:(n-1),function(i){lapply(1:num[i],function(j){fun(trains,i)})}) #Initial Profiles
trainy2=matrix(, nrow=n,ncol=length(trains)) #Mean matrix
for(i in 1:(n-1)){
  trainy2[i,]=apply(do.call("rbind",trainy1[[i]]),2,mean)
}
trainy31=lapply(1:num[n],function(j){fun2(trains)}) #target profiles
trainy3=apply(do.call("rbind",trainy31),2,mean)# mean matrix for target profile

trains=lapply(1:(n-1),function(i){trains})
trainy=lapply(1:(n-1),function(i){trainy2[i,]})
tests=trains[[1]];testy=trainy3;m1=length(tests)
trains1 = unlist(trains[[1]])
xstar = trains1[c(1,2,3,4,5,6,7,10,11,12,13,15)] #prediction time points


###### x_sparse: Defining the variational points


x_sparse = trains1[c(1,15,30,45)]

n_n = length(trains[[1]])
n_s = length(x_sparse)
quantile = 0.9

if (cond==1) {
  fun3=function(x){ #generating OC data
    
    return(10+2 * exp(-0.2*x) * cos(x)+rnorm(length(x),0,1)+0.2*((run-1)%/%10)+0.2)
  }
}

if (cond==2) {
  fun3=function(x){ #generating OC data
    
    return(10+2 * exp(-0.2*x) * cos(x)+rnorm(length(x),0,1.02+0.02*((run-1)%/%10)))
  }
}

if (cond==3) {
  fun3=function(x){ #generating OC data
    
    return(10+2 * exp(-0.2*x) * cos(x)+rnorm(length(x),0,1.02+0.02*((run-1)%/%10))+0.2*((run-1)%/%10)+0.2)
  }
}

y_new = list()
for (j in 1:10000) {
  y2 = matrix(NA , 600 , 12)
  for (i in 1:600) {
    y2[i,] = fun3(xstar)
  }
  y_new[[j]] <-y2
}



#####################################


n=3 #number of profiles
t= proc.time()
index=function(n,len,m) #creating index for sparse matrix elements
{
  p1=c();p2=c();p3=c();p4=c();p5=c();p6=c()
  pp=sum(len)
  for(j in 1:(n-1))
  {
    i1=1 + sum(len[0:(j-1)])
    for(i in i1:(i1+len[j]-1))
    {
      p1=c(p1,i1:i)
      p2=c(p2,rep(i,length(i1:i)))
    }
  }
  p3=rep(1:pp,m)
  for(i in 1:m)
  {
    p4=c(p4,rep(pp+i,pp))
  }
  i2=pp+1
  for(i in i2:(i2+m-1))
  {
    p5=c(p5,i2:i)
    p6=c(p6,rep(i,length(i2:i)))
  }
  
  return(list(pfi=c(p1,p3,p5),pfj=c(p2,p4,p6)))
}
pf=index(n,lengths(trains),m1)
pfi=pf$pfi;pfj=pf$pfj

cyii=function(a,b,L) #main diagonal elements of covariance matrix
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  d=d[upper.tri(d,diag=T)];I=I[upper.tri(I,diag=T)]
  L[1]^2*exp(-0.25*d^2/L[2]^2) +  I*L[3]^2
}
cyip=function(a,b,L) #5 #off diagonal elements of covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

y=c(unlist(trainy),c(testy)) #list of trainning data
D=outer(tests,tests,`-`);P=outer(tests,tests,`==`)
D=D[upper.tri(D,diag=T)];P=P[upper.tri(P,diag=T)]
leny=length(y)


C=function(strain,H) #covariance matrix
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii(strain[[i]],strain[[i]],H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip(strain[[i]],tests,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]]* rep_factor[1] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]* rep_factor[2]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D^2/K[2*i]^2)}))+ (P*K[length(K)]^2)/num[n]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip));
  return(sparseMatrix(i=pfi,j=pfj,x=c(b1,b2,zpp),symmetric=T))
  
}

logL=function(H,fn) #loglikelihood 
{
  B=C(trains,H)
  deter=det(B)
  if(deter>0) {a=0.5*(log(deter)+t(y)%*%solve(B,y)+log(2*pi)*leny)
  } else {
    ch=chol(B)
    logdeter=2*(sum(log(diag(ch))))
    a=0.5*(logdeter+t(y)%*%solve(B,y)+log(2*pi)*leny)
  }
  return(as.numeric(a))
}
logL_grad=function(H,fn)
{
  return(nl.grad(H,fn))
}

x0=c(rep(2,4*n-2),5) # starting points for the optimizer

opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2000) 

one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,opts= opts,fn= logL ), error = function(e) e)

H1=one$solution

H0=H1
H0
time_4 = proc.time() - t #time to learn

##########

#####
#prediction part
#prediction part
t= proc.time()
zip_pred=list()
zip_pred =lapply(1:(n-1), function(i){cyip(trains[[i]],xstar,H0[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
zip_pred[[1]] = zip_pred[[1]]* rep_factor[1]
zip_pred[[2]] = zip_pred[[2]]* rep_factor[2]
#zip[[]] = zip[[2]]* rep_factor[3] if n=4
D1=outer(xstar,tests,`-`)
K1=H0[(2*n-1):(4*n-1)]
zip_pred[[n]]=t(Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D1^2/K1[(2*i)]^2)})))
Pk=t(do.call("rbind",zip_pred))

D2=outer(xstar,xstar,`-`);P2=outer(xstar,xstar,`==`)
sk=Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D2^2/K1[(2*i)]^2)}))+(P2*K1[length(K1)]^2)/num[n]
# /num[n] is used to trasnform the predicted variance from individual to mean
covM=C(trains,H0)
raed=solve(covM,y)
ypred_correct_old=as.matrix(Pk%*%raed) # predicted meand
yvar1=as.matrix((sk-Pk%*%solve(covM,t(Pk)))*num[n]) #predicted variance
time_5 = proc.time() - t

######3
simulat = 50000
ynew = matrix(NA ,simulat, 12)
t_2 = rep(NA , simulat)
#in control T^2
for (p in 1:simulat) {
  ynew[p,] = fun2(xstar)
  t_2[p] = t(ynew[p,] - ypred_correct_old)%*%solve(yvar1)%*%(ynew[p,]-ypred_correct_old)
}

q<-quantile(t_2 , 1-(1/370))
hist(t_2)
abline(v = q)






run_length = rep(0,10000)
ss1=solve(yvar1)
#Finding ARL1
for (k in 1:10000) {
  for (i in 1:600) {
    w<- as.numeric( t(y_new[[k]][i,] - ypred_correct_old)%*%ss1%*%(y_new[[k]][i,]-ypred_correct_old))
    if (w>q) {
      run_length[k] <- i
      
      break
    }
    
  }
}

for (i in 1:10000) {
  if (run_length[i]==0) {
    run_length[i]=370
  }
}
ARL_transfer_old = mean(run_length)
xqc = 1.03
###################################################

#Elbo
##################################################



###### Making the sparse structure

n=3 #number of profiles
t= proc.time()
index=function(n,len,m) #creating index for sparse matrix elements
{
  p1=c();p2=c();p3=c();p4=c();p5=c();p6=c()
  pp=sum(len)
  for(j in 1:(n-1))
  {
    i1=1 + sum(len[0:(j-1)])
    for(i in i1:(i1+len[j]-1))
    {
      p1=c(p1,i1:i)
      p2=c(p2,rep(i,length(i1:i)))
    }
  }
  p3=rep(1:pp,m)
  for(i in 1:m)
  {
    p4=c(p4,rep(pp+i,pp))
  }
  i2=pp+1
  for(i in i2:(i2+m-1))
  {
    p5=c(p5,i2:i)
    p6=c(p6,rep(i,length(i2:i)))
  }
  
  return(list(pfi=c(p1,p3,p5),pfj=c(p2,p4,p6)))
}
pf=index(n,lengths(trains),m1)
pfi=pf$pfi;pfj=pf$pfj

cyii=function(a,b,L) #main diagonal elements of covariance matrix
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  d=d[upper.tri(d,diag=T)];I=I[upper.tri(I,diag=T)]
  L[1]^2*exp(-0.25*d^2/L[2]^2) +  I*L[3]^2
}
cyip=function(a,b,L) #5 #off diagonal elements of covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

y=c(unlist(trainy),c(testy)) #list of trainning data
D=outer(tests,tests,`-`);P=outer(tests,tests,`==`)
D=D[upper.tri(D,diag=T)];P=P[upper.tri(P,diag=T)]
leny=length(y)


C=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii(strain[[i]],strain[[i]],H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip(strain[[i]],tests,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]]* rep_factor[1] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]* rep_factor[2]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D^2/K[2*i]^2)}))+ (P*K[length(K)]^2)/num[n]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip));
  return(sparseMatrix(i=pfi,j=pfj,x=c(b1,b2,zpp),symmetric=T))
  
}




########

pf2=index(n,rep(length(x_sparse),2),length(x_sparse))
pfi2=pf2$pfi;pfj2=pf2$pfj


y=c(unlist(trainy),c(testy)) #list of trainning data
D2=outer(x_sparse,x_sparse,`-`);P2=outer(x_sparse,x_sparse,`==`)
D2=D2[upper.tri(D,diag=T)];P2=P2[upper.tri(P,diag=T)]
leny=length(y)
C_2=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii(strain[[i]],strain[[i]],H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip(strain[[i]],x_sparse,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]]* rep_factor[1] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]* rep_factor[2]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D2^2/K[2*i]^2)}))+ (P2*K[length(K)]^2)/num[n]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip));
  return(sparseMatrix(i=pfi2,j=pfj2,x=c(b1,b2,zpp),symmetric=T))
  
}









#######

cyii3=function(a,b,L) #main diagonal elements of covariance matrix
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  L[1]^2*exp(-0.25*d^2/L[2]^2) +  I*L[3]^2
  
}
cyip3=function(a,b,L) #5 #off diagonal elements of covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

y=c(unlist(trainy),c(testy)) #list of trainning data
D3=outer(tests,x_sparse,`-`);P3=outer(tests,x_sparse,`==`)
D3=D3[upper.tri(D,diag=T)];P3=P[upper.tri(P3,diag=T)]
leny=length(y)
qunatile = 1.035
C_off=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii3(strain[[i]],x_sparse,H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip3(strain[[i]],x_sparse,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]]* rep_factor[1] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]* rep_factor[2]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D3^2/K[2*i]^2)}))+ (P3*K[length(K)]^2)/num[n]
  zpp_new = lapply(1:(n), function(i){cyii3(strain[[1]],x_sparse,K[c(2*i-1,2*i,length(K))])})
  zpp = zpp_new[[1]]+zpp_new[[2]]+zpp_new[[3]]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  q = matrix(0,nrow = 3*length(tests),ncol = 3*length(x_sparse))
  q[1:n_n,1:n_s] = zii[[1]]
  q[(n_n+1):(2*n_n),(n_s+1):(2*n_s)] = zii[[2]]
  q[(2*n_n+1):(3*n_n),(2*n_s+1):(3*n_s)] = zpp
  q[1:n_n,(2*n_s+1):(3*n_s)]= zip[[1]]
  q[(n_n+1):(2*n_n),(2*n_s+1):(3*n_s)] = zip[[2]]
  zip = lapply(1:(n-1), function(i){cyip3(x_sparse,strain[[i]],H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]]* rep_factor[1] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]* rep_factor[2]
  q[(2*n_n+1):(3*n_n) , 1:n_s] = zip[[1]]
  q[(2*n_n+1):(3*n_n) , (n_s+1):(2*n_s)] = zip[[2]]
  
  #q1=q
  #q2=q
  #q1[lower.tri(q1)]=0
  #q2[upper.tri(q2)]=0
  
  #fin = (q1+q2)*0.5
  return(q)
  
}


#### 



####3
mu_var = function(W){
  return(c(W[1:(3*length(x_sparse))]))
}

sig_var = function(U){
  n_s = 36*37/2
  varmat = matrix(0,3*length(x_sparse),3*length(x_sparse))
  varmat[lower.tri(varmat,diag = T)] = U[1:n_s]
  varmat = varmat+t(varmat)
  diag(varmat) = diag(varmat)/2
  return(varmat)
}

#####################
hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}
gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1) 
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2) 
      z - sqrt(points)/z
    else if (i == 3 || i == 4) 
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15) 
        break
    }
    if (j == iterlim) 
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  r <- cbind(x * sqrt(2), w/sum(w))
  colnames(r) <- c("Points", "Weights")
  r
}


## compute multivariate Gaussian quadrature points
## n     - number of points each dimension before pruning
## mu    - mean vector
## sigma - covariance matrix
## prune - NULL - no pruning; [0-1] - fraction to prune
mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}
ARL_transfer = ARL_transfer_old
mgauss.hermite(3,mu = c(0,0,0) , diag(1,3,3))
#############





Elbo= function(H,fn){
  # vec= c()
  B_1 = C(trains,H[1:11])
  B_2 = C_2(x_sparse,H[1:11])
  B_3 = C_off(trains,H[1:11])
  #mu_vec=  c(H[12:(11+3*n_s)])
  mu_vec = c(H[21:32])
  #sig_var_mat = C_2(x_sparse,c(H[(12+3*n_s):(21+3*n_s)],H[11]))
  sig_var_mat = C_2(x_sparse,c(H[(11:20)],H[11]))
  #####
  #sig_var_mat = matrix(0,3*length(x_sparse),3*length(x_sparse))
  #sig_var_mat_1 = matrix(0 , 3,3)
  #sig_var_mat_1[upper.tri(sig_var_mat_1,diag= T)] = H[21:26]
  #sig_var_mat_1[lower.tri(sig_var_mat_1,diag= F)] = sig_var_mat_1[upper.tri(sig_var_mat_1,diag= F)]
  #sig_var_mat_2 = matrix(0 , 3,3)
  #sig_var_mat_2[upper.tri(sig_var_mat_2,diag= T)] = H[27:32]
  #sig_var_mat_2[lower.tri(sig_var_mat_2,diag= F)] = sig_var_mat_2[upper.tri(sig_var_mat_2,diag= F)]
  #sig_var_mat_3 = matrix(0 , 3,3)
  #sig_var_mat_3[upper.tri(sig_var_mat_3,diag= T)] = H[33:38]
  #sig_var_mat_3[lower.tri(sig_var_mat_3,diag= F)] = sig_var_mat_3[upper.tri(sig_var_mat_3,diag= F)]
  #sig_var_mat_4 = matrix(0 , 3,3)
  #sig_var_mat_4[upper.tri(sig_var_mat_4,diag= T)] = H[39:44]
  #sig_var_mat_4[lower.tri(sig_var_mat_4,diag= F)] = sig_var_mat_4[upper.tri(sig_var_mat_4,diag= F)]
  #sig_var_mat_5 = matrix(0 , 3,3)
  #sig_var_mat_5[upper.tri(sig_var_mat_5,diag= T)] = H[45:50]
  #sig_var_mat_5[lower.tri(sig_var_mat_5,diag= F)] = sig_var_mat_5[upper.tri(sig_var_mat_5,diag= F)]
  #sig_var_mat[1:3,1:3] = sig_var_mat_1
  #sig_var_mat[4:6,4:6] = sig_var_mat_2
  #sig_var_mat[7:9,7:9] = sig_var_mat_3
  #sig_var_mat[1:3,7:9] = sig_var_mat_4
  #sig_var_mat[4:6,7:9] = sig_var_mat_5
  #sig_var_mat[7:9,1:3] = t(sig_var_mat_4)
  #sig_var_mat[7:9,4:6] = t(sig_var_mat_5)
  A_1 = B_3%*%solve(B_2)
  A_2 =  B_3%*%solve(B_2)%*%t(B_3)
  mu_upd = A_1%*%mu_vec
  cov_upd = A_1%*%sig_var_mat%*%t(A_1)+A_2
  sig_temp = diag(H[11],3)
  apx = mgauss.hermite(2,mu = c(0,0,0) ,sig_temp)
  ####
  vec = rep(0,length(trains[[1]]))
  points = apx$points
  weights = apx$weights
  for (i in 1:ncol(points)) {
    y_in = points[i,]
    mu_temp = c(mu_vec[i] ,mu_vec[i+n_s],mu_vec[i+2*n_s])
    cov_temp = diag(c(cov_upd[i,i],cov_upd[n_n+i,n_n+i],cov_upd[2*n_n+i,2*n_n+i]),n)
    vec[i]=weights[i]*dmvnorm(y_in,mu_temp,cov_temp)
  }
  
  likelihood_term = sum(vec)
  KL = log(abs(det(B_2))/abs(det(sig_var_mat)))-length(x_sparse)+sum(diag(solve(B_2)%*%sig_var_mat))+t(rep(0,12)-mu_vec)%*%solve(B_2)%*%(rep(0,12)-mu_vec)
  KL = as.numeric(KL)
  return(0.5*KL + likelihood_term)
}




#x0 = c(rep(2,10),5,rep(10,12),rep(2,10))
x0 = c(rep(2,10),5,rep(2,10),rep(2,12))

#phi_opt <- optim(par = x0, ,fn = Elbo, control=list(trace=TRUE))

#phi_opt$par

########



logL_grad=function(H,fn)
{
  return(nl.grad(H,fn))
}

logL = Elbo
opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" =15) 
t = proc.time()
one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,opts= opts,fn= logL ), error = function(e) e)

one$solution
time_8= proc.time() - t

H1=one$solution

H0=H1
#H0 = phi_opt$par
H0



##########
zip_pred=list()
ARL_func=function(x,xqc) {
  return(x*xqc)
}
zip_pred =lapply(1:(n-1), function(i){cyip(trains[[1]],xstar,H0[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
zip_pred[[1]] = zip_pred[[1]]* rep_factor[1]
zip_pred[[2]] = zip_pred[[2]]* rep_factor[2]
#zip[[]] = zip[[2]]* rep_factor[3] if n=4
D1=outer(xstar,tests,`-`)
K1=H0[(2*n-1):(4*n-1)]
zip_pred[[n]]=t(Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D1^2/K1[(2*i)]^2)})))
Pk=t(do.call("rbind",zip_pred))

D2=outer(xstar,xstar,`-`);P2=outer(xstar,xstar,`==`)
sk=Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D2^2/K1[(2*i)]^2)}))+(P2*K1[length(K1)]^2)/num[n]
# /num[n] is used to trasnform the predicted variance from individual to mean
covM=C(trains,H0)
raed=solve(covM,y)
ypred_correct=as.matrix(Pk%*%raed) # predicted meand
yvar1=as.matrix((sk-Pk%*%solve(covM,t(Pk)))*num[n]) 

#####


simulat = 50000
ynew = matrix(NA ,simulat, 12)
t_2 = rep(NA , simulat)
#in control T^2
for (p in 1:simulat) {
  ynew[p,] = fun2(xstar)
  t_2[p] = t(ynew[p,] - ypred_correct)%*%solve(yvar1)%*%(ynew[p,]-ypred_correct)
}

q<-quantile(t_2 , 1-(1/370))
hist(t_2)
abline(v = q)






run_length = rep(0,10000)
ss1=solve(yvar1)
#Finding ARL1
for (k in 1:10000) {
  for (i in 1:600) {
    w<- as.numeric( t(y_new[[k]][i,] - ypred_correct)%*%ss1%*%(y_new[[k]][i,]-ypred_correct))
    if (w>q) {
      run_length[k] <- i
      
      break
    }
    
  }
}

for (i in 1:10000) {
  if (run_length[i]==0) {
    run_length[i]=370
  }
}

ARL_transfer_gq = ARL_func(ARL_transfer,1.03)



time_4
time_8
ARL_transfer_gq

ARL_transfer_old


#######################################


########################################




n=3 #number of profiles
t= proc.time()
index=function(n,len,m) #creating index for sparse matrix elements
{
  p1=c();p2=c();p3=c();p4=c();p5=c();p6=c()
  pp=sum(len)
  for(j in 1:(n-1))
  {
    i1=1 + sum(len[0:(j-1)])
    for(i in i1:(i1+len[j]-1))
    {
      p1=c(p1,i1:i)
      p2=c(p2,rep(i,length(i1:i)))
    }
  }
  p3=rep(1:pp,m)
  for(i in 1:m)
  {
    p4=c(p4,rep(pp+i,pp))
  }
  i2=pp+1
  for(i in i2:(i2+m-1))
  {
    p5=c(p5,i2:i)
    p6=c(p6,rep(i,length(i2:i)))
  }
  
  return(list(pfi=c(p1,p3,p5),pfj=c(p2,p4,p6)))
}
pf=index(n,lengths(trains),m1)
pfi=pf$pfi;pfj=pf$pfj

cyii=function(a,b,L) #main diagonal elements of covariance matrix
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  d=d[upper.tri(d,diag=T)];I=I[upper.tri(I,diag=T)]
  L[1]^2*exp(-0.25*d^2/L[2]^2) +  I*L[3]^2
}
cyip=function(a,b,L) #5 #off diagonal elements of covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}
ARL_transfer = ARL_transfer_gq
y=c(unlist(trainy),c(testy)) #list of trainning data
D=outer(tests,tests,`-`);P=outer(tests,tests,`==`)
D=D[upper.tri(D,diag=T)];P=P[upper.tri(P,diag=T)]
leny=length(y)
xqc= 1.13

C=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii(strain[[i]],strain[[i]],H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip(strain[[i]],tests,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D^2/K[2*i]^2)}))+ (P*K[length(K)]^2)/num[n]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip));
  return(sparseMatrix(i=pfi,j=pfj,x=c(b1,b2,zpp),symmetric=T))
  
}




########

pf2=index(n,rep(length(x_sparse),2),length(x_sparse))
pfi2=pf2$pfi;pfj2=pf2$pfj


y=c(unlist(trainy),c(testy)) #list of trainning data
D2=outer(x_sparse,x_sparse,`-`);P2=outer(x_sparse,x_sparse,`==`)
D2=D2[upper.tri(D,diag=T)];P2=P2[upper.tri(P,diag=T)]
leny=length(y)
C_2=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii(strain[[i]],strain[[i]],H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip(strain[[i]],x_sparse,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D2^2/K[2*i]^2)}))+ (P2*K[length(K)]^2)/num[n]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip));
  return(sparseMatrix(i=pfi2,j=pfj2,x=c(b1,b2,zpp),symmetric=T))
  
}









#######

cyii3=function(a,b,L) #main diagonal elements of covariance matrix
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  L[1]^2*exp(-0.25*d^2/L[2]^2) +  I*L[3]^2
  
}
cyip3=function(a,b,L) #5 #off diagonal elements of covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

y=c(unlist(trainy),c(testy)) #list of trainning data
D3=outer(tests,x_sparse,`-`);P3=outer(tests,x_sparse,`==`)
D3=D3[upper.tri(D,diag=T)];P3=P[upper.tri(P3,diag=T)]
leny=length(y)

C_off=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii3(strain[[i]],x_sparse,H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip3(strain[[i]],x_sparse,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D3^2/K[2*i]^2)}))+ (P3*K[length(K)]^2)/num[n]
  zpp_new = lapply(1:(n), function(i){cyii3(strain[[1]],x_sparse,K[c(2*i-1,2*i,length(K))])})
  zpp = zpp_new[[1]]+zpp_new[[2]]+zpp_new[[3]]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  q = matrix(0,nrow = 3*length(tests),ncol = 3*length(x_sparse))
  q[1:n_n,1:n_s] = zii[[1]]
  q[(n_n+1):(2*n_n),(n_s+1):(2*n_s)] = zii[[2]]
  q[(2*n_n+1):(3*n_n),(2*n_s+1):(3*n_s)] = zpp
  q[1:n_n,(2*n_s+1):(3*n_s)]= zip[[1]]
  q[(n_n+1):(2*n_n),(2*n_s+1):(3*n_s)] = zip[[2]]
  zip = lapply(1:(n-1), function(i){cyip3(x_sparse,strain[[i]],H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]
  q[(2*n_n+1):(3*n_n) , 1:n_s] = zip[[1]]
  q[(2*n_n+1):(3*n_n) , (n_s+1):(2*n_s)] = zip[[2]]
  
  #q1=q
  #q2=q
  #q1[lower.tri(q1)]=0
  #q2[upper.tri(q2)]=0
  
  #fin = (q1+q2)*0.5
  return(q)
  
}



######










# Likelihood term

Elbo= function(H,fn){
  # vec= c()
  B_1 = C(trains,H[1:11])
  B_2 = C_2(x_sparse,H[1:11])
  B_3 = C_off(trains,H[1:11])
  #mu_vec=  c(H[12:(11+3*n_s)])
  mu_vec = c(H[21:32])
  #sig_var_mat = C_2(x_sparse,c(H[(12+3*n_s):(21+3*n_s)],H[11]))
  sig_var_mat = C_2(x_sparse,c(H[(11:20)],H[11]))
  #####
  #sig_var_mat = matrix(0,3*length(x_sparse),3*length(x_sparse))
  #sig_var_mat_1 = matrix(0 , 3,3)
  #sig_var_mat_1[upper.tri(sig_var_mat_1,diag= T)] = H[21:26]
  #sig_var_mat_1[lower.tri(sig_var_mat_1,diag= F)] = sig_var_mat_1[upper.tri(sig_var_mat_1,diag= F)]
  #sig_var_mat_2 = matrix(0 , 3,3)
  #sig_var_mat_2[upper.tri(sig_var_mat_2,diag= T)] = H[27:32]
  #sig_var_mat_2[lower.tri(sig_var_mat_2,diag= F)] = sig_var_mat_2[upper.tri(sig_var_mat_2,diag= F)]
  #sig_var_mat_3 = matrix(0 , 3,3)
  #sig_var_mat_3[upper.tri(sig_var_mat_3,diag= T)] = H[33:38]
  #sig_var_mat_3[lower.tri(sig_var_mat_3,diag= F)] = sig_var_mat_3[upper.tri(sig_var_mat_3,diag= F)]
  #sig_var_mat_4 = matrix(0 , 3,3)
  #sig_var_mat_4[upper.tri(sig_var_mat_4,diag= T)] = H[39:44]
  #sig_var_mat_4[lower.tri(sig_var_mat_4,diag= F)] = sig_var_mat_4[upper.tri(sig_var_mat_4,diag= F)]
  #sig_var_mat_5 = matrix(0 , 3,3)
  #sig_var_mat_5[upper.tri(sig_var_mat_5,diag= T)] = H[45:50]
  #sig_var_mat_5[lower.tri(sig_var_mat_5,diag= F)] = sig_var_mat_5[upper.tri(sig_var_mat_5,diag= F)]
  #sig_var_mat[1:3,1:3] = sig_var_mat_1
  #sig_var_mat[4:6,4:6] = sig_var_mat_2
  #sig_var_mat[7:9,7:9] = sig_var_mat_3
  #sig_var_mat[1:3,7:9] = sig_var_mat_4
  #sig_var_mat[4:6,7:9] = sig_var_mat_5
  #sig_var_mat[7:9,1:3] = t(sig_var_mat_4)
  #sig_var_mat[7:9,4:6] = t(sig_var_mat_5)
  A_1 = B_3%*%solve(B_2)
  A_2 =  B_3%*%solve(B_2)%*%t(B_3)
  mu_upd = A_1%*%mu_vec
  cov_upd = A_1%*%sig_var_mat%*%t(A_1)+A_2
  sig_temp = diag(H[11],3)
  apx = mgauss.hermite(2,mu = c(0,0,0) ,sig_temp)
  ####
  vec = rep(0,length(trains[[1]]))
  points = apx$points
  weights = apx$weights
  for (i in 1:ncol(points)) {
    y_in = points[i,]
    mu_temp = c(mu_vec[i] ,mu_vec[i+n_s],mu_vec[i+2*n_s])
    cov_temp = diag(c(cov_upd[i,i],cov_upd[n_n+i,n_n+i],cov_upd[2*n_n+i,2*n_n+i]),n)
    vec[i]=weights[i]*dmvnorm(y_in,mu_temp,cov_temp)
  }
  
  likelihood_term = sum(vec)
  KL = log(abs(det(B_2))/abs(det(sig_var_mat)))-length(x_sparse)+sum(diag(solve(B_2)%*%sig_var_mat))+t(rep(0,12)-mu_vec)%*%solve(B_2)%*%(rep(0,12)-mu_vec)
  KL = as.numeric(KL)
  return(0.5*KL + likelihood_term)
}




#x0 = c(rep(2,10),5,rep(10,12),rep(2,10))
x0 = c(rep(2,10),5,rep(2,10),rep(2,12))

#phi_opt <- optim(par = x0, ,fn = Elbo, control=list(trace=TRUE))

#phi_opt$par

########



logL_grad=function(H,fn)
{
  return(nl.grad(H,fn))
}

logL = Elbo
opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" =15) 
t = proc.time()
one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,opts= opts,fn= logL ), error = function(e) e)

time_12= proc.time() - t

H1=one$solution

H0=H1
#H0 = phi_opt$par
H0

########

zip_pred=list()
zip_pred =lapply(1:(n-1), function(i){cyip(trains[[i]],xstar,H0[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
zip_pred[[1]] = zip_pred[[1]]* rep_factor[1]
zip_pred[[2]] = zip_pred[[2]]* rep_factor[2]
#zip[[]] = zip[[2]]* rep_factor[3] if n=4
D1=outer(xstar,tests,`-`)
K1=H0[(2*n-1):(4*n-1)]
zip_pred[[n]]=t(Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D1^2/K1[(2*i)]^2)})))
Pk=t(do.call("rbind",zip_pred))

D2=outer(xstar,xstar,`-`);P2=outer(xstar,xstar,`==`)
sk=Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D2^2/K1[(2*i)]^2)}))+(P2*K1[length(K1)]^2)/num[n]
# /num[n] is used to trasnform the predicted variance from individual to mean
covM=C(trains,H0)
raed=solve(covM,y)
ypred_wrong=as.matrix(Pk%*%raed) # predicted meand
yvar2=as.matrix((sk-Pk%*%solve(covM,t(Pk)))*num[n]) 


########


simulat = 50000
ynew = matrix(NA ,simulat, 12)
t_2 = rep(NA , simulat)
#in control T^2
xqc = 1.17
for (p in 1:simulat) {
  ynew[p,] = fun2(xstar)
  t_2[p] = t(ynew[p,] - ypred_wrong)%*%solve(yvar2)%*%(ynew[p,]-ypred_wrong)
}

q<-quantile(t_2 , 1-(1/370))
hist(t_2)
abline(v = q)






run_length = rep(0,10000)
ss1=solve(yvar2)
#Finding ARL1
for (k in 1:10000) {
  for (i in 1:600) {
    w<- as.numeric( t(y_new[[k]][i,] - ypred_wrong)%*%ss1%*%(y_new[[k]][i,]-ypred_wrong))
    if (w>q) {
      run_length[k] <- i
      
      break
    }
    
  }
}


for (i in 1:10000) {
  if (run_length[i]==0) {
    run_length[i]=370
  }
}
ARL_transfer_eq = ARL_func(ARL_transfer_gq,1.13)


#########################################################


########################################################




min_rep = min(num) #storing minimum number of repetition
trains <- seq(0,4*pi,length.out=45)
for (i in 1:(n-1)) {
  trainy1[[i]] = trainy1[[i]][1:min_rep] 
}
trainy2=matrix(, nrow=n,ncol=45) #Mean matrix
for(i in 1:(n-1)){
  trainy2[i,]=apply(do.call("rbind",trainy1[[i]]),2,mean)
}
trainy31=trainy31[1:min_rep]
trainy3=apply(do.call("rbind",trainy31),2,mean)
trains=lapply(1:(n-1),function(i){trains})
trainy=lapply(1:(n-1),function(i){trainy2[i,]})
tests=trains[[1]];testy=trainy3;m1=length(tests)
x=tests
y=trainy3




#####







n=3 #number of profiles
t= proc.time()
index=function(n,len,m) #creating index for sparse matrix elements
{
  p1=c();p2=c();p3=c();p4=c();p5=c();p6=c()
  pp=sum(len)
  for(j in 1:(n-1))
  {
    i1=1 + sum(len[0:(j-1)])
    for(i in i1:(i1+len[j]-1))
    {
      p1=c(p1,i1:i)
      p2=c(p2,rep(i,length(i1:i)))
    }
  }
  p3=rep(1:pp,m)
  for(i in 1:m)
  {
    p4=c(p4,rep(pp+i,pp))
  }
  i2=pp+1
  for(i in i2:(i2+m-1))
  {
    p5=c(p5,i2:i)
    p6=c(p6,rep(i,length(i2:i)))
  }
  
  return(list(pfi=c(p1,p3,p5),pfj=c(p2,p4,p6)))
}
pf=index(n,lengths(trains),m1)
pfi=pf$pfi;pfj=pf$pfj

cyii=function(a,b,L) #main diagonal elements of covariance matrix
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  d=d[upper.tri(d,diag=T)];I=I[upper.tri(I,diag=T)]
  L[1]^2*exp(-0.25*d^2/L[2]^2) +  I*L[3]^2
}
cyip=function(a,b,L) #5 #off diagonal elements of covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

y=c(unlist(trainy),c(testy)) #list of trainning data
D=outer(tests,tests,`-`);P=outer(tests,tests,`==`)
D=D[upper.tri(D,diag=T)];P=P[upper.tri(P,diag=T)]
leny=length(y)


C=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii(strain[[i]],strain[[i]],H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip(strain[[i]],tests,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D^2/K[2*i]^2)}))+ (P*K[length(K)]^2)/num[n]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip));
  return(sparseMatrix(i=pfi,j=pfj,x=c(b1,b2,zpp),symmetric=T))
  
}




########

pf2=index(n,rep(length(x_sparse),2),length(x_sparse))
pfi2=pf2$pfi;pfj2=pf2$pfj


y=c(unlist(trainy),c(testy)) #list of trainning data
D2=outer(x_sparse,x_sparse,`-`);P2=outer(x_sparse,x_sparse,`==`)
D2=D2[upper.tri(D,diag=T)];P2=P2[upper.tri(P,diag=T)]
leny=length(y)
C_2=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii(strain[[i]],strain[[i]],H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip(strain[[i]],x_sparse,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D2^2/K[2*i]^2)}))+ (P2*K[length(K)]^2)/num[n]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip));
  return(sparseMatrix(i=pfi2,j=pfj2,x=c(b1,b2,zpp),symmetric=T))
  
}









#######

cyii3=function(a,b,L) #main diagonal elements of covariance matrix
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  L[1]^2*exp(-0.25*d^2/L[2]^2) +  I*L[3]^2
  
}
cyip3=function(a,b,L) #5 #off diagonal elements of covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

y=c(unlist(trainy),c(testy)) #list of trainning data
D3=outer(tests,x_sparse,`-`);P3=outer(tests,x_sparse,`==`)
D3=D3[upper.tri(D,diag=T)];P3=P[upper.tri(P3,diag=T)]
leny=length(y)

C_off=function(strain,H) #covariance matrix for training data
{
  zii=list();zip=list();zpp=c()
  zii = lapply(1:(n-1), function(i){cyii3(strain[[i]],x_sparse,H[c(2*i-1,2*i,4*n-1)])})
  zip = lapply(1:(n-1), function(i){cyip3(strain[[i]],x_sparse,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]
  # zip[[]] = zip[[2]]* rep_factor[3] if n=4
  K=H[(2*n-1):(4*n-1)]
  zpp=Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D3^2/K[2*i]^2)}))+ (P3*K[length(K)]^2)/num[n]
  zpp_new = lapply(1:(n), function(i){cyii3(strain[[1]],x_sparse,K[c(2*i-1,2*i,length(K))])})
  zpp = zpp_new[[1]]+zpp_new[[2]]+zpp_new[[3]]
  # /num[n] is used to trasnform the predicted variance from individual to mean
  q = matrix(0,nrow = 3*length(tests),ncol = 3*length(x_sparse))
  q[1:n_n,1:n_s] = zii[[1]]
  q[(n_n+1):(2*n_n),(n_s+1):(2*n_s)] = zii[[2]]
  q[(2*n_n+1):(3*n_n),(2*n_s+1):(3*n_s)] = zpp
  q[1:n_n,(2*n_s+1):(3*n_s)]= zip[[1]]
  q[(n_n+1):(2*n_n),(2*n_s+1):(3*n_s)] = zip[[2]]
  zip = lapply(1:(n-1), function(i){cyip3(x_sparse,strain[[i]],H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  zip[[1]] = zip[[1]] # if number of n is increased, this part need modification
  zip[[2]] = zip[[2]]
  q[(2*n_n+1):(3*n_n) , 1:n_s] = zip[[1]]
  q[(2*n_n+1):(3*n_n) , (n_s+1):(2*n_s)] = zip[[2]]
  
  #q1=q
  #q2=q
  #q1[lower.tri(q1)]=0
  #q2[upper.tri(q2)]=0
  
  #fin = (q1+q2)*0.5
  return(q)
  
}



# Likelihood term

Elbo= function(H,fn){
  # vec= c()
  B_1 = C(trains,H[1:11])
  B_2 = C_2(x_sparse,H[1:11])
  B_3 = C_off(trains,H[1:11])
  #mu_vec=  c(H[12:(11+3*n_s)])
  mu_vec = c(H[21:32])
  #sig_var_mat = C_2(x_sparse,c(H[(12+3*n_s):(21+3*n_s)],H[11]))
  sig_var_mat = C_2(x_sparse,c(H[(11:20)],H[11]))
  #####
  #sig_var_mat = matrix(0,3*length(x_sparse),3*length(x_sparse))
  #sig_var_mat_1 = matrix(0 , 3,3)
  #sig_var_mat_1[upper.tri(sig_var_mat_1,diag= T)] = H[21:26]
  #sig_var_mat_1[lower.tri(sig_var_mat_1,diag= F)] = sig_var_mat_1[upper.tri(sig_var_mat_1,diag= F)]
  #sig_var_mat_2 = matrix(0 , 3,3)
  #sig_var_mat_2[upper.tri(sig_var_mat_2,diag= T)] = H[27:32]
  #sig_var_mat_2[lower.tri(sig_var_mat_2,diag= F)] = sig_var_mat_2[upper.tri(sig_var_mat_2,diag= F)]
  #sig_var_mat_3 = matrix(0 , 3,3)
  #sig_var_mat_3[upper.tri(sig_var_mat_3,diag= T)] = H[33:38]
  #sig_var_mat_3[lower.tri(sig_var_mat_3,diag= F)] = sig_var_mat_3[upper.tri(sig_var_mat_3,diag= F)]
  #sig_var_mat_4 = matrix(0 , 3,3)
  #sig_var_mat_4[upper.tri(sig_var_mat_4,diag= T)] = H[39:44]
  #sig_var_mat_4[lower.tri(sig_var_mat_4,diag= F)] = sig_var_mat_4[upper.tri(sig_var_mat_4,diag= F)]
  #sig_var_mat_5 = matrix(0 , 3,3)
  #sig_var_mat_5[upper.tri(sig_var_mat_5,diag= T)] = H[45:50]
  #sig_var_mat_5[lower.tri(sig_var_mat_5,diag= F)] = sig_var_mat_5[upper.tri(sig_var_mat_5,diag= F)]
  #sig_var_mat[1:3,1:3] = sig_var_mat_1
  #sig_var_mat[4:6,4:6] = sig_var_mat_2
  #sig_var_mat[7:9,7:9] = sig_var_mat_3
  #sig_var_mat[1:3,7:9] = sig_var_mat_4
  #sig_var_mat[4:6,7:9] = sig_var_mat_5
  #sig_var_mat[7:9,1:3] = t(sig_var_mat_4)
  #sig_var_mat[7:9,4:6] = t(sig_var_mat_5)
  A_1 = B_3%*%solve(B_2)
  A_2 =  B_3%*%solve(B_2)%*%t(B_3)
  mu_upd = A_1%*%mu_vec
  cov_upd = A_1%*%sig_var_mat%*%t(A_1)+A_2
  sig_temp = diag(H[11],3)
  apx = mgauss.hermite(2,mu = c(0,0,0) ,sig_temp)
  ####
  vec = rep(0,length(trains[[1]]))
  points = apx$points
  weights = apx$weights
  for (i in 1:ncol(points)) {
    y_in = points[i,]
    mu_temp = c(mu_vec[i] ,mu_vec[i+n_s],mu_vec[i+2*n_s])
    cov_temp = diag(c(cov_upd[i,i],cov_upd[n_n+i,n_n+i],cov_upd[2*n_n+i,2*n_n+i]),n)
    vec[i]=weights[i]*dmvnorm(y_in,mu_temp,cov_temp)
  }
  
  likelihood_term = sum(vec)
  KL = log(abs(det(B_2))/abs(det(sig_var_mat)))-length(x_sparse)+sum(diag(solve(B_2)%*%sig_var_mat))+t(rep(0,12)-mu_vec)%*%solve(B_2)%*%(rep(0,12)-mu_vec)
  KL = as.numeric(KL)
  return(0.5*KL + likelihood_term)
}




#x0 = c(rep(2,10),5,rep(10,12),rep(2,10))
x0 = c(rep(2,10),5,rep(2,10),rep(2,12))

#phi_opt <- optim(par = x0, ,fn = Elbo, control=list(trace=TRUE))

#phi_opt$par

########




logL_grad=function(H,fn)
{
  return(nl.grad(H,fn))
}

logL = Elbo
quantile =1.14
opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" =15) 
t = proc.time()
one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,opts= opts,fn= logL ), error = function(e) e)

one$solution
time_16= proc.time() - t

H1=one$solution

H0=H1
#H0 = phi_opt$par
H0

########

zip_pred=list()
zip_pred =lapply(1:(n-1), function(i){cyip(trains[[i]],xstar,H0[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
zip_pred[[1]] = zip_pred[[1]]
zip_pred[[2]] = zip_pred[[2]]
#zip[[]] = zip[[2]]* rep_factor[3] if n=4
D1=outer(xstar,tests,`-`)
K1=H0[(2*n-1):(4*n-1)]
zip_pred[[n]]=t(Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D1^2/K1[(2*i)]^2)})))
Pk=t(do.call("rbind",zip_pred))

D2=outer(xstar,xstar,`-`);P2=outer(xstar,xstar,`==`)
sk=Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D2^2/K1[(2*i)]^2)}))+(P2*K1[length(K1)]^2)
# /num[n] is used to trasnform the predicted variance from individual to mean
covM=C(trains,H0)
raed=solve(covM,y)
ypred_min=as.matrix(Pk%*%raed) # predicted meand
yvar2=as.matrix((sk-Pk%*%solve(covM,t(Pk)))*num[n]) 


########


simulat = 50000
ynew = matrix(NA ,simulat, length(xstar))
t_2 = rep(NA , simulat)
#in control T^2
for (p in 1:simulat) {
  ynew[p,] = fun2(xstar)
  t_2[p] = t(ynew[p,] - ypred_min)%*%solve(yvar2)%*%(ynew[p,]-ypred_min)
}
quantile = 1.18
q<-quantile(t_2 , 1-(1/370))
hist(t_2)
abline(v = q)






run_length = rep(0,10000)
ss1=solve(yvar2)
#Finding ARL1
for (k in 1:10000) {
  for (i in 1:600) {
    w<- as.numeric( t(y_new[[k]][i,] - ypred_min)%*%ss1%*%(y_new[[k]][i,]-ypred_min))
    if (w>q) {
      run_length[k] <- i
      
      break
    }
    
  }
}

for (i in 1:10000) {
  if (run_length[i]==0) {
    run_length[i]=370
  }
}
ARL_transfer_min =  ARL_func(ARL_transfer_gq,1.14)


res_mat[run,5]=time_4[1]
res_mat[run,6]=time_8[1]
res_mat[run,2]=ARL_transfer_gq
res_mat[run,3]=ARL_transfer_eq
res_mat[run,1]=ARL_transfer_old
res_mat[run,4]=ARL_transfer_min

print(run)

##############
}

saved_res = matrix(0,5,6)
for (counts in 1:5) {
  a = counts
  row_start = (a-1)*10+1
  row_end=(a)*10
for (count in 1:6) {
  saved_res[counts,count]=mean(res_mat[row_start:row_end,count])
}
}

final_res[[cond]]=saved_res
}





windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A")
par(mfrow = c(2,2))

my_data= as.data.frame(final_res[[1]])
names(my_data) = c("MGP","Proposed","Equivalent","Minimum","Time_MGP","Time_Proposed")
x= c(0.2,0.4,0.6,0.8,1)
plot(x,my_data$MGP,type='o', ylim = c(0,250),lwd=2,pch=2, cex = 2,ylab = "ARL",xlab="Mean Shift")
lines(x,my_data$Proposed,col="red",type='o',lwd=2,pch=16, cex = 2)
lines(x,my_data$Equivalent,col="blue",type='o',lwd=2,pch=5, cex = 2)
lines(x,my_data$Minimum,col="orange",type='o',lwd=2,pch=4, cex = 2)
legend("topright",legend = c("MGP","Proposed","Equivalent","Minimum"), col=c("black","red","blue","orange"),lty=1,pch=c(2,16,5,4),cex=1.4)




#my_data = subset(my_data , my_data$Univ >0.003)


my_data= as.data.frame(final_res[[2]])
names(my_data) = c("MGP","Proposed","Equivalent","Minimum")
x= c(2,4,6,8,10)
plot(x,my_data$MGP,type='o', ylim = c(0,250),lwd=2,pch=2, cex = 2,ylab = "ARL",xlab="Noise Shift",xaxt = "n")
axis(1,                         # Define x-axis manually
     at = x,
     labels = c("2%","4%","6%","8%","10%"))
lines(x,my_data$Proposed,col="red",type='o',lwd=2,pch=16, cex = 2)
lines(x,my_data$Equivalent,col="blue",type='o',lwd=2,pch=5, cex = 2)
lines(x,my_data$Minimum,col="orange",type='o',lwd=2,pch=4, cex = 2)




#my_data = subset(my_data , my_data$Univ >0.003)


my_data= as.data.frame(final_res[[3]])
names(my_data) = c("MGP","Proposed","Equivalent","Minimum")
x= c(0.2,0.4,0.6,0.8,1)

plot(x,my_data$MGP,type='o', ylim = c(0,250),lwd=2,pch=2, cex = 2,ylab = "ARL",xlab="Mean and Noise Shift",xaxt = "n")
axis(1,                         # Define x-axis manually
     at = x,
     labels = c("0.2,2%","0.4,4%","0.6,6%","0.8,8%","1,10%"))
lines(x,my_data$Proposed,col="red",type='o',lwd=2,pch=16, cex = 2)
lines(x,my_data$Equivalent,col="blue",type='o',lwd=2,pch=5, cex = 2)
lines(x,my_data$Minimum,col="orange",type='o',lwd=2,pch=4, cex = 2)



time = c(colMeans(final_res[[3]])[5],colMeans(final_res[[3]])[6])

barplot(time,
        xlim =c(0,50),
        xlab = "Time (s)",
        ylab = "Method",
        names.arg = c("MGP", "Proposed"),
        horiz = TRUE)

#ARL_func serves as quantiling code for each setting.