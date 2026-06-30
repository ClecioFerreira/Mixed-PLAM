#########################################################################################################################
#########################################################################################################################
######################## AR Skew Mixed
#################### Residuo quantilico tq_i

tqiSMSN<-function(object){ 
  data <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  depStruct <- object$depStruct
  timeVar <- object$timeVar
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  p<-ncol(x)
  q1<-ncol(z)
  N<- nrow(data)
  ind_levels <- levels(ind)
  distr <- object$distr
  
  beta1 <- as.matrix(object$estimates$beta)
  sigmae <- as.numeric(object$estimates$sigma2)
  lambda=as.matrix(object$estimates$lambda)
  D <-object$estimates$D
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Delta = matrix.sqrt(D)%*%delta
  xsi=matrix.sqrt(solve(D))%*%lambda

  alphas <- object$alphas
  nknots <- object$nknots
  Nm <- object$Nmatrix
  qnl <- length(Nm)
  gammas <- object$estimates$gammas
  
  if (distr=="sn"|distr=="norm") {c.=-sqrt(2/pi)}
  if (distr=="st"|distr=="t") {c.=-sqrt(object$estimates$nu/pi)*
    gamma((object$estimates$nu-1)/2)/gamma(object$estimates$nu/2)}
  if (distr=="ssl"|distr=="sl") {c.=-sqrt(2/pi)*object$estimates$nu/(object$estimates$nu-.5)}
  if (distr=="scn"|distr=="cn") {c.=-sqrt(2/pi)*(1+object$estimates$nu[1]*(object$estimates$nu[2]^(-.5)-1))
  }
  
  tqi=rep(0,length(ind_levels))
  index <- 1:N
  for (i in seq_along(ind_levels)) {
    seqi <- index[ind==ind_levels[i]]
    timei <- time[ind==ind_levels[i]]
    xi <- matrix(x[seqi,],ncol=p)
    zi <- matrix(z[seqi,],ncol=q1)
    yi <- as.matrix(y[seqi])
    lower = rep(-Inf,length(seqi))
    if (depStruct=="UNC") Sigmai <- sigmae*diag(length(seqi))
    if (depStruct=="ARp") Sigmai <- sigmae*CovARp(object$estimates$phi,timei)
    if (depStruct=="CS") Sigmai <- sigmae*CovCS(object$estimates$phi,length(seqi))
    if (depStruct=="DEC") Sigmai <- sigmae*CovDEC(object$estimates$phi[1],object$estimates$phi[2],timei)
    if (depStruct=="CAR1") Sigmai <- sigmae*CovDEC(object$estimates$phi,1,timei)
    medg=0
    for (j in 1:qnl){
      Ni=Nm[[j]][seqi,]
      mui=Ni%*%gammas[[j]]
      medg=medg+mui
    }
    mui = xi%*%beta1 +medg + c.*zi%*%Delta
    Psii = Sigmai + zi%*%D%*%t(zi)
    Lambdai=solve(solve(D)+t(zi)%*%solve(Sigmai)%*%zi)
    lambdai=matrix.sqrt(solve(Psii))%*%zi%*%D%*%xsi/sqrt(1+as.numeric(t(xsi)%*%Lambdai%*%xsi))
    if (distr=="sn") acum <- pmvSN(lower,yi,mui,Psii,lambdai)
    if (distr=="st") acum <- pmvST(lower,yi,mui,Psii,lambdai,nu=object$estimates$nu)
    if (distr=="scn") {
      nu=object$estimates$nu
      acum <- nu[1]*pmvSN(lower,yi,mui,Psii/nu[2],lambdai)+(1-nu[1])*pmvSN(lower,yi,mui,Psii,lambdai)
    }
    tqi[i]=qnorm(acum)
  }
  return(tqi)
}


tqiSMN<-function(object){ 
  data <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  depStruct <- object$depStruct
  timeVar <- object$timeVar
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  p<-ncol(x)
  q1<-ncol(z)
  N<- nrow(data)
  ind_levels <- levels(ind)
  distr <- object$distr
  
  beta1 <- as.matrix(object$estimates$beta)
  sigmae <- as.numeric(object$estimates$sigma2)
  lambda=matrix(0,q1,1)
  D <-object$estimates$D
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Delta = matrix.sqrt(D)%*%delta
  Gamma = D - Delta%*%t(Delta)
  xsi=matrix.sqrt(solve(D))%*%lambda
  
  alphas <- object$alphas
  nknots <- object$nknots
  Nm <- object$Nmatrix
  qnl <- length(Nm)
  gammas <- object$estimates$gammas
  
  if (distr=="sn"|distr=="norm") {c.=-sqrt(2/pi)}
  if (distr=="st"|distr=="t") {c.=-sqrt(object$estimates$nu/pi)*
    gamma((object$estimates$nu-1)/2)/gamma(object$estimates$nu/2)}
  if (distr=="ssl"|distr=="sl") {c.=-sqrt(2/pi)*object$estimates$nu/(object$estimates$nu-.5)}
  if (distr=="scn"|distr=="cn") {c.=-sqrt(2/pi)*(1+object$estimates$nu[1]*(object$estimates$nu[2]^(-.5)-1))
  }
  
  tqi=rep(0,length(ind_levels))
  index <- 1:N
  for (i in seq_along(ind_levels)) {
    seqi <- index[ind==ind_levels[i]]
    timei <- time[ind==ind_levels[i]]
    xi <- matrix(x[seqi,],ncol=p)
    zi <- matrix(z[seqi,],ncol=q1)
    yi <- as.matrix(y[seqi])
    lower = rep(-Inf,length(seqi))
    if (depStruct=="UNC") Sigmai <- sigmae*diag(length(seqi))
    if (depStruct=="ARp") Sigmai <- sigmae*CovARp(object$estimates$phi,timei)
    if (depStruct=="CS") Sigmai <- sigmae*CovCS(object$estimates$phi,length(seqi))
    if (depStruct=="DEC") Sigmai <- sigmae*CovDEC(object$estimates$phi[1],object$estimates$phi[2],timei)
    if (depStruct=="CAR1") Sigmai <- sigmae*CovDEC(object$estimates$phi,1,timei)
    medg=0
    for (j in 1:qnl){
      Ni=Nm[[j]][seqi,]
      mui=Ni%*%gammas[[j]]
      medg=medg+mui
    }
    mui = xi%*%beta1 +medg + c.*zi%*%Delta
    Psii = Sigmai + zi%*%D%*%t(zi)
    Lambdai=solve(solve(D)+t(zi)%*%solve(Sigmai)%*%zi)
    lambdai=matrix.sqrt(solve(Psii))%*%zi%*%D%*%xsi/sqrt(1+as.numeric(t(xsi)%*%Lambdai%*%xsi))
    if (distr=="norm") acum <- pmvSN(lower,yi,mui,Psii,lambdai)
    if (distr=="t") acum <- pmvST(lower,yi,mui,Psii,lambdai,nu=object$estimates$nu)
    if (distr=="cn") {
      nu=object$estimates$nu
      acum <- nu[1]*pmvSN(lower,yi,mui,Psii/nu[2],lambdai)+(1-nu[1])*pmvSN(lower,yi,mui,Psii,lambdai)
    }
    tqi[i]=qnorm(acum)
  }
  return(tqi)
}


envel<-function(tqi,alpha=0.05){
  replic <- 1000
  n <- length(tqi)
  aa<-sort(tqi,index.return=TRUE)
  tqi2<-aa$x
  ordem=aa$ix
  e <- matrix(0,n,replic)
  e1 <- numeric(n)
  e2 <- numeric(n)
  #
  for(i in 1:replic) e[,i] <- sort(rnorm(n,0,1))
  #
  for(i in 1:n){
    eo <- sort(e[i,])
    e1[i] <- eo[ceiling(alpha/2*replic)]
    e2[i] <- eo[ceiling((1-alpha/2)*replic)]
    
    }
  #
  med <- apply(e,1,mean)
  faixa <- range(tqi,e1,e2)
  #
  aux<-0
  for(i in 1:n){
    if ((tqi2[i]<e1[i])|(tqi2[i]>e2[i])) aux<-c(aux,i) # |=ou
  }
  #print(aux)
  ll=length(aux)
  aux2=aux[2:ll]
  if (ll >= 1) posi=ordem[aux2]
  #
  par(pty="s")
  qqnorm(tqi,xlab="Theoretical standard normal quantiles",
         ylab="Sample values and simulated envelope", ylim=faixa, pch=16,font.lab=2)
  par(new=T)
  qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
  par(new=T)
  qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
  par(new=T)
  qqnorm(med,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=2)
  #--------------------------------------------------------------#
  envel2<-posi
  print(length(posi))
  print(posi)
  
}



resMahalSMSN<-function(object){ 
  data <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  depStruct <- object$depStruct
  timeVar <- object$timeVar
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  p<-ncol(x)
  q1<-ncol(z)
  N<- nrow(data)
  ind_levels <- levels(ind)
  distr <- object$distr
  
  beta1 <- as.matrix(object$estimates$beta)
  sigmae <- as.numeric(object$estimates$sigma2)
  lambda=matrix(0,q1,1)
  if (distr=="sn"|distr=="st"|distr=="ssl"|distr=="scn") {lambda=as.matrix(object$estimates$lambda)}
  D <-object$estimates$D
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Delta = matrix.sqrt(D)%*%delta
  xsi=matrix.sqrt(solve(D))%*%lambda
  
  alphas <- object$alphas
  nknots <- object$nknots
  Nm <- object$Nmatrix
  qnl <- length(Nm)
  gammas <- object$estimates$gammas
  
  if (distr=="sn"|distr=="norm") {c.=-sqrt(2/pi)}
  if (distr=="st"|distr=="t") {c.=-sqrt(object$estimates$nu/pi)*
    gamma((object$estimates$nu-1)/2)/gamma(object$estimates$nu/2)}
  if (distr=="ssl"|distr=="sl") {c.=-sqrt(2/pi)*object$estimates$nu/(object$estimates$nu-.5)}
  if (distr=="scn"|distr=="cn") {c.=-sqrt(2/pi)*(1+object$estimates$nu[1]*(object$estimates$nu[2]^(-.5)-1))
  }
  
  resi=rep(0,length(ind_levels))
  index <- 1:N
  for (i in seq_along(ind_levels)) {
    seqi <- index[ind==ind_levels[i]]
    timei <- time[ind==ind_levels[i]]
    xi <- matrix(x[seqi,],ncol=p)
    zi <- matrix(z[seqi,],ncol=q1)
    yi <- as.matrix(y[seqi])
    lower = rep(-Inf,length(seqi))
    if (depStruct=="UNC") Sigmai <- sigmae*diag(length(seqi))
    if (depStruct=="ARp") Sigmai <- sigmae*CovARp(object$estimates$phi,timei)
    if (depStruct=="CS") Sigmai <- sigmae*CovCS(object$estimates$phi,length(seqi))
    if (depStruct=="DEC") Sigmai <- sigmae*CovDEC(object$estimates$phi[1],object$estimates$phi[2],timei)
    if (depStruct=="CAR1") Sigmai <- sigmae*CovDEC(object$estimates$phi,1,timei)
    medg=0
    for (j in 1:qnl){
      Ni=Nm[[j]][seqi,]
      mui=Ni%*%gammas[[j]]
      medg=medg+mui
    }
    mui = xi%*%beta1 +medg + c.*zi%*%Delta
    Psii = Sigmai + zi%*%D%*%t(zi)
    #Lambdai=solve(solve(D)+t(zi)%*%solve(Sigmai)%*%zi)
    #lambdai=matrix.sqrt(solve(Psii))%*%zi%*%D%*%xsi/sqrt(1+as.numeric(t(xsi)%*%Lambdai%*%xsi))
    Mahal_i=as.numeric(t(yi-mui)%*%solve(Psii)%*%(yi-mui))
    ni=length(seqi)
    if (distr=="sn"|distr=="norm") {
      mui1=1-2/(9*ni)
      vari=2/(9*ni)
      T=(Mahal_i/ni)^(1/3)
      resi[i]=(T-mui1)/sqrt(vari)
      }
    if (distr=="st"|distr=="t") {
      Fi=Mahal_i/ni
      nu=object$estimates$nu
      resi[i]=((1-2/(9*nu))*(Fi^(1/3))-(1-2/(9*ni)))/sqrt(2/(9*nu)*(Fi^(2/3))+2/(9*ni))
    }
    if (distr=="scn"|distr=="cn") {
      nu=object$estimates$nu
      nu1=nu[1]
      nu2=nu[2]
      Fi=Mahal_i
      mu=ni*nu1/(nu2^2)+ni*(1-nu1)
      var=2*ni*(nu1^2)/(nu2^4)+2*ni*(1-nu1)^2
      resi[i]=(Fi-mu)/sqrt(var)
    }
    
  }
  return(resi)
}



envel_boot<-function(object,alpha=0.05){
  replic <- 1000
  data <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  depStruct <- object$depStruct
  timeVar <- object$timeVar
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  p<-ncol(x)
  q1<-ncol(z)
  N<- nrow(data)
  ind_levels <- levels(ind)
  distr <- object$distr
  
  beta1 <- as.matrix(object$estimates$beta)
  sigmae <- as.numeric(object$estimates$sigma2)
  lambda=matrix(0,q1,1)
  if (distr=="sn"|distr=="st"|distr=="ssl"|distr=="scn") {lambda=as.matrix(object$estimates$lambda)}
  D <-object$estimates$D
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Delta = matrix.sqrt(D)%*%delta
  xsi=matrix.sqrt(solve(D))%*%lambda
  phi=object$estimates$phi
  nu=NULL
  if (distr=="t"|distr=="st"|distr=="sl"|distr=="ssl"|distr=="cn"|distr=="scn") {nu=object$estimates$nu}
  Nm <- object$Nmatrix
  qnl <- length(Nm)
  gammas <- object$estimates$gammas
  if (distr=="sn"|distr=="norm") {c.=-sqrt(2/pi)}
  if (distr=="st"|distr=="t") {c.=-sqrt(object$estimates$nu/pi)*
    gamma((object$estimates$nu-1)/2)/gamma(object$estimates$nu/2)}
  if (distr=="ssl"|distr=="sl") {c.=-sqrt(2/pi)*object$estimates$nu/(object$estimates$nu-.5)}
  if (distr=="scn"|distr=="cn") {c.=-sqrt(2/pi)*(1+object$estimates$nu[1]*(object$estimates$nu[2]^(-.5)-1))
  }
  n=length(ind_levels)
  Xsim <- c()# matrix(0,n,replic)
  index <- 1:N
  Mahal_obs=rep(0,n)
  for (i in seq_along(ind_levels)) {
    seqi <- index[ind==ind_levels[i]]
    timei <- time[ind==ind_levels[i]]
    xi <- matrix(x[seqi,],ncol=p)
    zi <- matrix(z[seqi,],ncol=q1)
    yi <- as.matrix(y[seqi])
    
    if (depStruct=="UNC") Sigmai <- sigmae*diag(length(seqi))
    if (depStruct=="ARp") Sigmai <- sigmae*CovARp(object$estimates$phi,timei)
    if (depStruct=="CS") Sigmai <- sigmae*CovCS(object$estimates$phi,length(seqi))
    if (depStruct=="DEC") Sigmai <- sigmae*CovDEC(object$estimates$phi[1],object$estimates$phi[2],timei)
    if (depStruct=="CAR1") Sigmai <- sigmae*CovDEC(object$estimates$phi,1,timei)
    
    medg=0
    for (j in 1:qnl){
      Ni=Nm[[j]][seqi,]
      mui=Ni%*%gammas[[j]]
      medg=medg+mui
    }
    mui=xi%*%beta1+medg
    muic = xi%*%beta1 +medg + c.*zi%*%Delta
    Psii = Sigmai + zi%*%D%*%t(zi)
    Psi_inv=solve(Psii)
    Mahal_obs[i]=as.numeric(t(yi-muic)%*%Psi_inv%*%(yi-muic))
    Mahal=rep(0,replic)
    ni=length(seqi)
    for (k in 1:replic){
    auxi=rsmsn.lmm(1:ni,mui,zi,sigma2e,D,beta1,lambda,depStruct=depStruct,phi=phi,distr=distr,nu=nu)
    Yk=matrix(auxi$y)
    Mahal[k]=as.numeric(t(Yk-muic)%*%Psi_inv%*%(Yk-muic))
    }
    Xsim <- cbind(Xsim,Mahal)
  }
  ###################
  d2med <-apply(Xsim,2,mean)
  d21<-rep(0,n)
  d22<-rep(0,n)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim[,i],alpha/2)
    d22[i]  <- quantile(Xsim[,i],1-alpha/2)
  }
  ###########
  aa<-sort(Mahal_obs,index.return=TRUE)
  Mahal2<-aa$x
  ordem=aa$ix
  aux<-0
  for(i in 1:n){
    if ((Mahal2[i]<d21[i])|(Mahal2[i]>d22[i])) aux<-c(aux,i) # |=ou
  }
  ll=length(aux)
  aux2=aux[2:ll]
  if (ll >= 1) posi=ordem[aux2]

  ##########################
  print(ordem)
print(d21[ordem])
  print(d2med[ordem])
  
  fy <- range(Mahal2,d21,d22)
  xq2 <- qchisq(ppoints(n), 1)
  plot(xq2,Mahal2,xlab = expression(bold(paste("Theoretical ",chi[1]^2, " quantile"))),
       ylab="Sample value and simulated envelope",pch=20,ylim=fy,font.lab=2)
  #par(new=T)
  lines(xq2,d21[ordem])#,type="l",ylim=fy,xlab="",ylab="")
  #par(new=T)
  lines(xq2,d2med[ordem])#,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
  #par(new=T)
  lines(xq2,d22[ordem])#,type="l",ylim=fy,xlab="",ylab="")
  print(posi)
}



plot_log_density<-function(object,replic=1000){# OHagan et al 2017
  data <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  depStruct <- object$depStruct
  timeVar <- object$timeVar
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  p<-ncol(x)
  q1<-ncol(z)
  N<- nrow(data)
  ind_levels <- levels(ind)
  distr <- object$distr
  alphas <- object$alphas
  Km <- object$Kmatrix
  
  beta1 <- as.matrix(object$estimates$beta)
  sigmae <- as.numeric(object$estimates$sigma2)
  lambda=matrix(0,q1,1)
  if (distr=="sn"|distr=="st"|distr=="ssl"|distr=="scn") {lambda=as.matrix(object$estimates$lambda)}
  D <-object$estimates$D
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Delta = matrix.sqrt(D)%*%delta
  Gamma=D-Delta%*%t(Delta)
  xsi=matrix.sqrt(solve(D))%*%lambda
  phi=object$estimates$phi
  nu=NULL
  if (distr=="t"|distr=="st"|distr=="sl"|distr=="ssl"|distr=="cn"|distr=="scn") {nu=object$estimates$nu}
  Nm <- object$Nmatrix
  qnl <- length(Nm)
  gammas <- object$estimates$gammas
  if (distr=="sn"|distr=="norm") {c.=-sqrt(2/pi)}
  if (distr=="st"|distr=="t") {c.=-sqrt(object$estimates$nu/pi)*
    gamma((object$estimates$nu-1)/2)/gamma(object$estimates$nu/2)}
  if (distr=="ssl"|distr=="sl") {c.=-sqrt(2/pi)*object$estimates$nu/(object$estimates$nu-.5)}
  if (distr=="scn"|distr=="cn") {c.=-sqrt(2/pi)*(1+object$estimates$nu[1]*(object$estimates$nu[2]^(-.5)-1))
  }
  n=length(ind_levels)
  Xsim <- c()# matrix(0,n,replic)
  index <- 1:N
  log_dens=rep(0,n)
  log_densb=matrix(0,n,replic)
  for (i in seq_along(ind_levels)) {
    seqi <- index[ind==ind_levels[i]]
    ni=length(seqi)
    timei <- time[ind==ind_levels[i]]
    xi <- matrix(x[seqi,],ncol=p)
    zi <- matrix(z[seqi,],ncol=q1)
    yi <- as.matrix(y[seqi])
    
    if (depStruct=="UNC") Sigmai <- sigmae*diag(length(seqi))
    if (depStruct=="ARp") Sigmai <- sigmae*CovARp(object$estimates$phi,timei)
    if (depStruct=="CS") Sigmai <- sigmae*CovCS(object$estimates$phi,length(seqi))
    if (depStruct=="DEC") Sigmai <- sigmae*CovDEC(object$estimates$phi[1],object$estimates$phi[2],timei)
    if (depStruct=="CAR1") Sigmai <- sigmae*CovDEC(object$estimates$phi,1,timei)
    
    medg=0
    penalty=0
    for (j in 1:qnl){
      gammai=matrix(gammas[[j]])
      Ni=Nm[[j]][seqi,]
      mui=Ni%*%gammai
      medg=medg+mui
      penaltyi=alphas[j]/ni*t(gammai)%*%Km[[j]]%*%gammai
      penalty=penalty+as.numeric(penaltyi)
    }
    mui=xi%*%beta1+medg
    muic = xi%*%beta1 +medg + c.*zi%*%Delta
    resi=yi-muic
    Psii = Sigmai + zi%*%D%*%t(zi)
    Psi_inv=solve(Psii)
    Psim=solve(matrix.sqrt(Psii))
    di<-as.numeric(t(resi)%*%Psi_inv%*%(resi))
    Mtj2 = (1+t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%zi%*%Delta)^(-1)
    Ajj  = as.numeric(sqrt(Mtj2)*t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%resi)
    
    if (distr=="norm") log_dens[i]=log(dmvnorm(as.numeric(yi),as.numeric(muic),Psii)) -0.5/n*penalty
    if (distr=="t"){
      dti = gamma((nu+ni)/2)/gamma(nu/2)/pi^(ni/2)/sqrt(det(Psii))*nu^(-ni/2)*(di/nu+1)^(-(ni+nu)/2)
      log_dens[i]=log(dti) -0.5/n*penalty
    }
    if (distr=="sl"){
      f2 <- function(u) u^(nu - 1)*((2*pi)^(-ni/2))*(u^(ni/2))*((det(Psii))^(-1/2))*exp(-0.5*u*di)
      resp <- integrate(Vectorize(f2),0,1)$value
      log_dens[i]=log(nu*resp) - 0.5/n*penalty
    }
    if (distr=="cn"){
      log_dens[i]= log(nu[1]*dmvnorm(as.numeric(yi),as.numeric(muic),(Psii/nu[2]))+ (1-nu[1])*dmvnorm(as.numeric(yi),as.numeric(muic),Psii)) - 0.5/n*penalty
    }
    if (distr=="sn"){
      log_dens[i]=log(2*dmvnorm(as.numeric(yi),as.numeric(muic),Psii)*pnorm(Ajj,0,1))-0.5/n*penalty
    }
    if (distr=="st"){
      #Mtj2 = (1+t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%zi%*%Delta)^(-1)
      #Ajj = as.numeric(sqrt(Mtj2)*t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%resi)
      dtj = gamma((nu+ni)/2)/gamma(nu/2)/pi^(ni/2)/sqrt(det(Psii))*nu^(-ni/2)*(di/nu+1)^(-(ni+nu)/2)
      log_dens[i]=log(2*dtj*pt(sqrt(nu+ni)*Ajj/sqrt(di+nu),nu+ni))-0.5/n*penalty
    }
    if (distr=="ssl"){
      #Mtj2 = (1+t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%zi%*%Delta)^(-1)
      #Ajj  = as.numeric(sqrt(Mtj2)*t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%resi)
      f2 = function(u) u^(nu - 1)*((2*pi)^(-ni/2))*(u^(ni/2))*((det(Psii))^(-1/2))*exp(-0.5*u*di)*pnorm(u^(1/2)*Ajj)
      resp = integrate(Vectorize(f2),0,1)$value
      log_dens[i]=log(2*nu*resp) -0.5/n*penalty
    }
    if (distr=="scn"){
      #Mtj2 = (1+t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%zi%*%Delta)^(-1)
      #Ajj  = as.numeric(sqrt(Mtj2)*t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%resi)
      log_dens[i]=log(2*(nu[1]*dmvnorm(as.numeric(yi),as.numeric(muic),(Psii/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
               (1-nu[1])*dmvnorm(as.numeric(yi),as.numeric(muic),Psii)*pnorm(Ajj,0,1)))  - 0.5/n*penalty
    }
    
    for (k in 1:replic){
      auxi=rsmsn.lmm(1:ni,mui,zi,sigma2e,D,beta1,lambda,depStruct=depStruct,phi=phi,distr=distr,nu=nu)
      yk=matrix(auxi$y)
      resi=yk-muic
      di<-as.numeric(t(resi)%*%Psi_inv%*%(resi))
      Mtj2 = (1+t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%zi%*%Delta)^(-1)
      Ajj  = as.numeric(sqrt(Mtj2)*t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%resi)
      
      if (distr=="norm") log_densb[i,k]=log(dmvnorm(as.numeric(yk),as.numeric(muic),Psii)) -0.5/n*penalty
      if (distr=="t"){
        dti = gamma((nu+ni)/2)/gamma(nu/2)/pi^(ni/2)/sqrt(det(Psii))*nu^(-ni/2)*(di/nu+1)^(-(ni+nu)/2)
        log_densb[i,k]=log(dti) -0.5/n*penalty
      }
      if (distr=="sl"){
        f2 <- function(u) u^(nu - 1)*((2*pi)^(-ni/2))*(u^(ni/2))*((det(Psii))^(-1/2))*exp(-0.5*u*di)
        resp <- integrate(Vectorize(f2),0,1)$value
        log_densb[i,k]=log(nu*resp) - 0.5/n*penalty
      }
      if (distr=="cn"){
        log_densb[i,k]= log(nu[1]*dmvnorm(as.numeric(yk),as.numeric(muic),(Psii/nu[2]))+ (1-nu[1])*dmvnorm(as.numeric(yk),as.numeric(muic),Psii)) - 0.5/n*penalty
      }
      if (distr=="sn"){
        log_densb[i,k]=log(2*dmvnorm(as.numeric(yk),as.numeric(muic),Psii)*pnorm(Ajj,0,1))-0.5/n*penalty
      }
      if (distr=="st"){
        #Mtj2 = (1+t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%zi%*%Delta)^(-1)
        #Ajj = as.numeric(sqrt(Mtj2)*t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%resi)
        dtj = gamma((nu+ni)/2)/gamma(nu/2)/pi^(ni/2)/sqrt(det(Psii))*nu^(-ni/2)*(di/nu+1)^(-(ni+nu)/2)
        log_densb[i,k]=log(2*dtj*pt(sqrt(nu+ni)*Ajj/sqrt(di+nu),nu+ni))-0.5/n*penalty
      }
      if (distr=="ssl"){
        #Mtj2 = (1+t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%zi%*%Delta)^(-1)
        #Ajj  = as.numeric(sqrt(Mtj2)*t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%resi)
        f2 = function(u) u^(nu - 1)*((2*pi)^(-ni/2))*(u^(ni/2))*((det(Psii))^(-1/2))*exp(-0.5*u*di)*pnorm(u^(1/2)*Ajj)
        resp = integrate(Vectorize(f2),0,1)$value
        log_densb[i,k]=log(2*nu*resp) -0.5/n*penalty
      }
      if (distr=="scn"){
        #Mtj2 = (1+t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%zi%*%Delta)^(-1)
        #Ajj  = as.numeric(sqrt(Mtj2)*t(Delta)%*%t(zi)%*%solve(Sigmai+zi%*%Gamma%*%t(zi))%*%resi)
        log_densb[i,k]=log(2*(nu[1]*dmvnorm(as.numeric(yk),as.numeric(muic),(Psii/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
                             (1-nu[1])*dmvnorm(as.numeric(yk),as.numeric(muic),Psii)*pnorm(Ajj,0,1)))  - 0.5/n*penalty
      }
    }
}
# Fim do for i e for k
  mtotal=cbind(log_dens,log_densb)
  xlims=c(min(mtotal),max(mtotal))

  plot(ecdf(log_densb[1,]),col="grey",xlab="Log density", ylab="ECDF",pch=16,cex=1,main="",xlim=xlims)
  par(new=T)
  for (k in 2:replic) {
    par(new=T)
    plot(ecdf(log_densb[,k]),col="grey",xlab="", ylab="",pch=16,cex=1,main="",yaxt="n", xaxt="n")}
  par(new=T)
  plot(ecdf(log_dens),xlab="", ylab="",pch=16,cex=1,main="",yaxt="n", xaxt="n")

  
}
