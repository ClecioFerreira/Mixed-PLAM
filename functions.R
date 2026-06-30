
require(nlme)
require(splines)
require(mvtnorm)
require(mnormt)
require(matrixcalc)
require(pracma)
require(MASS)
require(linpk)# pra usar blockdiag, mas soh funciona para submatrizes quadradas
require(Matrix) # para criar matriz bloco-diagonal as.matrix(Z))
require(future)
require(furrr)
require(tidytable)
require(MomTrunc)
#library(foreach)
#library(doParallel)
require(optimParallel)
source('auxfunctions.R')
source('auxfunctions-sim.R')
source('fixptfunctions.R')
source('objfunctions.R')
source('auxfunctions-daarem.R')
source('residuals.R')

skewness<-function(x){
   z=x-mean(x)
   skewness=mean(z^3)/((mean(z^2))^1.5)
   return(skewness)
}

kurtosis <- function (x, na.rm = FALSE) 
{
    var = function(x, ...) {
        mean((x - mean(x, ...))^2)
    }
    if (na.rm) 
        x <- x[!is.na(x)]
    sum((x - mean(x))^4)/(length(x) * var(x, na.rm = na.rm)^2) - 
        3
}

n_distinct<-function(ind) length(as.numeric(table(ind)))

# Construcao da matriz N
### Eilers and Marx
bsplinec<-function(x,K,bdeg=4){
# bdeg eh a ordem do polinomio; cubico default   
# d=bdeg-1, grau do polinomio
# K: numero de knots (interiores): 
# x: pontos para calcular matriz B do spline, sendo a=min(x) e b=max(x), limites dos dados; a=xsi_{d+1}; tau_1= xsi_{d+2}
  d=bdeg-1
  a=min(x)
  b=max(x)
  xsi1=seq(a,b,length=K+2)
  delta=xsi1[2]-xsi1[1]
  xsi=c(a-delta*(d:1),xsi1,b+delta*(1:d))
  B=splineDesign(xsi,x,bdeg,0*x) 
  return(B) 
}


is.wholenumber <- function(x, tol1 = .Machine$double.eps^0.5)  abs(x - round(x)) < tol1


matrix.sqrt <- function(A){
  if (length(A)==1) return(sqrt(A))
  else{
    sva <- svd(A)
    if (min(sva$d)>=0) {
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v) # svd e decomposi??o espectral
      if (all(abs(Asqrt%*%Asqrt-A)<1e-4)) return(Asqrt)
      else stop("Matrix square root is not defined/not real")
    }
    else stop("Matrix square root is not defined/not real")
  }
}

################################################################
# Trace of a matrix of dim >=1
################################################################
traceM <- function(Mat){
  if(length(Mat)==1) tr<- as.numeric(Mat)
  else tr<-sum(diag(Mat))
  tr
}

################################################################
# Inverter ordem de hierarquia de uma lista com nomes
################################################################
revert_list <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list)
}

# Transformation function: pi to phi
estphit <- function(pit) {
  p <- length(pit)
  Phi <- matrix(0,ncol=p,nrow=p)
  if (p>1) {
    diag(Phi) <- pit
    for (j in 2:p) {
      for (k in 1:(j-1)) {
        Phi[j,k] <- Phi[j-1,k] - pit[j]*Phi[j-1,j-k]
      }
    }
    return(Phi[p,])
  }
  else return(pit)
}

# Transformation function: phi to pi
tphitopi <- function(phit) {
  p <- length(phit)
  Phi <- matrix(0,ncol=p,nrow=p)
  Phi[p,] <- phit
  if (p>1) {
    for (k in p:2) {
      for (i in 1:(k-1)) {
        Phi[k-1,i] <- (Phi[k,i] + Phi[k,k]*Phi[k,k-i])/(1-Phi[k,k]^2)
      }
    }
    return(diag(Phi))
  }
  else return(phit)
}

# Matrix Mn = 1/sigma2 * autocov(arp) in function of vector of times
CovARp <- function(phi,ti) {
  p <- length(phi)
  n <- max(ti)
  if (n==1) Rn <- matrix(1)
  else Rn <- toeplitz(ARMAacf(ar=phi, ma=0, lag.max = n-1))
  rhos <- ARMAacf(ar=phi, ma=0, lag.max = p)[-1]
  Rn <- Rn/(1-sum(rhos*phi))
  return(as.matrix(Rn[ti,ti]))
}

# Corr matrix of comp symmetry (CS)
CovCS <- function(phi, n) {
  if (n==1) Rn <- matrix(1)
  else Rn <- toeplitz(c(1,rep(phi,n-1)))
  return(Rn)
}

# Corr matrix of CAR1
# CovDEC(phi,theta=1)

# Corr matrix of DEC
CovDEC <- function(phi1, phi2, ti) {
  ni <- length(ti)
  Rn <- diag(ni)
  if (ni==1) Rn <- matrix(1)
  else {
    for (i in 1:(ni-1)) for (j in (i+1):ni) Rn[i,j] <- phi1^(abs(ti[i]-ti[j])^phi2)
    Rn[lower.tri(Rn)] <- t(Rn)[lower.tri(Rn)]
  }
  return(Rn)
}

Dmatrix <- function(dd) {
  q2 <- length(dd)
  q1 <- -0.5+sqrt(1+8*q2)/2
  if (q1%%1 != 0) stop("wrong dimension of dd")
  D1 <- matrix(nrow=q1, ncol=q1)
  D1[upper.tri(D1,diag=T)] <- as.numeric(dd)
  D1[lower.tri(D1)] <- t(D1)[lower.tri(D1)]
  return(D1)
}

# Gerando smsn para um individuo usando rep hierarquica
gerar_ind_smsn = function(ni, Sig, Di, beta, lambda, distr="sn", nu=NULL) {
  if (distr=="sn") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st") {ui=rgamma(1,nu/2,nu/2)
                    c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {ui=rbeta(1,nu,1)
                    c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {ui=ifelse(runif(1)<nu[1],nu[2],1)
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = matrix.sqrt(Di)%*%delta
  Gammab = Di - Delta%*%t(Delta)
  Xi = cbind(1,runif(ni,0,2))
  Zi = matrix(1,nrow=ni)
  Beta = matrix(beta,ncol=1)
  q1 = nrow(Di)
  ti = c.+abs(rnorm(1,0,ui^-.5))
  bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,Xi%*%Beta+Zi%*%bi,sigma=ui^(-1)*Sig))

  return(data.frame(y=Yi,x=Xi[,2],tempo=1:ni,ui=ui))
}

gerar_ind_smsn_semip = function(ni, mui, Zi, Sig, D1, lambda, distr="sn", nu=NULL) {
# Geral: fazer usuario gerar mu_i
  if (distr=="sn") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st") {ui=rgamma(1,nu/2,nu/2)
                    c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {ui=rbeta(1,nu,1)
                    c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {ui=ifelse(runif(1)<nu[1],nu[2],1)
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = matrix.sqrt(D1)%*%delta
  Gammab = D1 - Delta%*%t(Delta)
  q1 = nrow(D1)
  ti = c.+abs(rnorm(1,0,ui^-.5))
  bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,mui+Zi%*%bi,sigma=ui^(-1)*Sig)) 
  return(data.frame(y=Yi,bi=bi,tempo=1:ni,ui=ui))
}

errorVar<- function(times,object=NULL,sigma2=NULL,depStruct=NULL,phi=NULL) {
  if((!is.null(object))&&(!inherits(object,c("SMSN","SMN")))) stop("object must inherit from class SMSN or SMN")
  if (is.null(object)&&is.null(depStruct)) stop("object or depStruct must be provided")
  if (is.null(object)&&is.null(sigma2)) stop("object or sigma2 must be provided")
  if (is.null(depStruct)) depStruct<-object$depStruct
  if (depStruct=="CI") depStruct = "UNC"
  if (depStruct!="UNC" && is.null(object)&is.null(phi)) stop("object or phi must be provided")
  if (!(depStruct %in% c("UNC","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or CAR1")
  if (is.null(sigma2)) sigma2<-object$estimates$sigma2
  if (is.null(phi)&&depStruct!="UNC") phi<-object$estimates$phi
  if (depStruct=="ARp" && (any(!is.wholenumber(times))|any(times<=0))) stop("times must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp" && any(tphitopi(phi)< -1|tphitopi(phi)>1)) stop("AR(p) non stationary, choose other phi")
  #
  if (depStruct=="UNC") var.out<- sigma2*diag(length(times))
  if (depStruct=="ARp") var.out<- sigma2*CovARp(phi,times)
  if (depStruct=="CS") var.out<- sigma2*CovCS(phi,length(times))
  if (depStruct=="DEC") var.out<- sigma2*CovDEC(phi[1],phi[2],times)
  if (depStruct=="CAR1") var.out<- sigma2*CovDEC(phi,1,times)
  var.out
}

rsmsn.lmm <- function(time1,mu,z1,sigma2,D1,beta,lambda,depStruct="UNC",phi=NULL,distr="sn",nu=NULL) {
  if (length(D1)==1 && !is.matrix(D1)) D1=as.matrix(D1)
  q1 = nrow(D1)
  p = length(beta)
  #if (ncol(as.matrix(x1))!=p) stop("incompatible dimension of x1/beta")
  if (ncol(as.matrix(z1))!=q1) stop ("incompatible dimension of z1/D1")
  if (length(lambda)!=q1) stop ("incompatible dimension of lambda/D1")
  if (!is.matrix(D1)) stop("D must be a matrix")
  if ((ncol(D1)!=q1)|(nrow(D1)!=q1)) stop ("wrong dimension of D")
  if (length(sigma2)!=1) stop ("wrong dimension of sigma2")
  if (sigma2<=0) stop("sigma2 must be positive")
  Sig <- errorVar(time1,depStruct = depStruct,sigma2=sigma2,phi=phi)
  #
  if (distr=="ssl") distr<-"ss"
  if (!(distr %in% c("norm","sn","t","st","sl","ss","cn","scn"))) stop("Invalid distribution")
  if (distr=="sn"|distr=="norm") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st"|distr=="t") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss"|distr=="sl") {ui=rbeta(1,nu,1); c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn"|distr=="cn") {ui=ifelse(runif(1)<nu[1],nu[2],1);
                      c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = matrix.sqrt(D1)%*%delta
  #Gammab = D1 - Delta%*%t(Delta)
  Zi = matrix(z1,ncol=q1)
  #ti = c.+abs(rnorm(1,0,ui^-.5))
  #bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  #Yi = t(rmvnorm(1,mu+Zi%*%bi,sigma=ui^(-1)*Sig))
  T0=abs(rnorm(1))
  T1=t(rmvnorm(1,rep(0,q1),diag(q1)))
  bi=c.*Delta+ui^(-0.5)*matrix.sqrt(D1)%*%(delta*T0+matrix.sqrt(diag(q1)-delta%*%t(delta))%*%T1 )
  Yi=mu+Zi%*%bi+t(rmvnorm(1,rep(0,ncol(Sig)),Sig))
  #if (all(Zi[,1]==1)) Zi = Zi[,-1]
  return(data.frame(time=time1,y=Yi))
}

lmmControl <- function(tol=1e-6,max.iter=300,calc.se=FALSE,
                       lb=NULL,lu=NULL,luDEC=3,
                       initialValues =list(beta=NULL,sigma2=NULL,D=NULL, lambda=NULL,phi=NULL,nu=NULL),
                       quiet=FALSE,showCriterium=FALSE,algorithm="DAAREM",
                       parallelphi=NULL, parallelnu=NULL, ncores=NULL,
                       control.daarem=list(),alphas=NULL,nknots=NULL) {
  if ((!is.numeric(tol))||(length(tol)>1)) stop("tol must be a small number")
  if (max.iter%%1!=0) {
    max.iter <- ceiling(max.iter)
    warning("using max.iter = ", max.iter)
  }
  if (!is.logical(calc.se)) stop("calc.se must be TRUE or FALSE")
  if ((!is.numeric(luDEC))||luDEC<0) stop("luDEC must be a (not too small) number")
  if (!is.list(initialValues)) stop("initialValues must be a list")
  if (all(c("D","dsqrt") %in% names(initialValues))) initialValues$dsqrt<-NULL
  if (any(!(names(initialValues) %in% c("beta","sigma2","lambda","D","phi","nu")))) warning("initialValues must be a list with named elements beta, sigma2, D, phi and/or nu, elements with other names are ignored")
  if (!is.list(control.daarem)) stop("control.daarem must be a list")
  if (!(algorithm%in%c("DAAREM","EM"))) stop("algorithm must be either 'EM' or 'DAAREM'")
  if (!is.logical(quiet)) stop("quiet must be TRUE or FALSE")
  if (!is.logical(showCriterium)) stop("showCriterium must be TRUE or FALSE")
  if (!is.null(parallelphi)) if (!is.logical(parallelphi)) stop("parallelphi must be TRUE or FALSE")
  if (!is.null(parallelnu)) if (!is.logical(parallelnu)) stop("parallelnu must be TRUE or FALSE")
  if (!is.null(ncores)) if (ncores%%1!=0) {
    ncores <- ceiling(ncores)
    warning("using ncores = ", ncores)
  }
  out<-list(tol = tol, max.iter = max.iter, calc.se = calc.se,
            lb = lb, lu = lu, luDEC = luDEC, initialValues = initialValues,
            quiet = quiet, showCriterium = showCriterium, algorithm = algorithm,
            parallelphi = parallelphi, parallelnu = parallelnu, ncores = ncores,
            control.daarem = control.daarem,alphas=alphas,nknots=nknots)
  class(out) <- c("lmmControl","list")
  return(out)
}


###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
############# SMSN models

smsn.lmm  <- function(data, formFixed, formFixedNL, groupVar, formRandom = ~1, depStruct = "UNC", timeVar = NULL, distr = "sn", covRandom = "pdSymm", skewind, pAR = 1, control = lmmControl())    # data tem de ser formato dataframe
{
    if (!is(formFixed, "formula")) 
        stop("formFixed must be a formula")
    if (!is(formRandom, "formula")) 
        stop("formRandom must be a formula")
    if (!inherits(control, "lmmControl")) 
        stop("control must be a list generated with lmmControl()")
    if (!is.character(groupVar)) 
        stop("groupVar must be a character containing the name of the grouping variable in data")
    if (!is.null(timeVar) && !is.character(timeVar)) 
        stop("timeVar must be a character containing the name of the time variable in data")
    if (length(formFixed) != 3) 
        stop("formFixed must be a two-sided linear formula object")
    if (!is.data.frame(data)) 
        stop("data must be a data.frame")
    if (length(class(data)) > 1) 
        data <- as.data.frame(data)
    vars_used <- unique(c(all.vars(formFixed), all.vars(formRandom), 
        groupVar, timeVar))
    vars_miss <- which(!(vars_used %in% names(data)))
    if (length(vars_miss) > 0) 
        stop(paste(vars_used[vars_miss], "not found in data"))
    data = data[, vars_used]
    if (!is.factor(data[, groupVar])) 
        data[, groupVar] <- haven::as_factor(data[, groupVar])
    data$ind <- data[, groupVar]
    ind <- data[, groupVar]
    x <- model.matrix(formFixed, data = data)
    y <- data[, all.vars(formFixed)[1]]
    z <- model.matrix(formRandom, data = data)
    m <- nlevels(ind)
    if (m <= 1) 
        stop(paste(groupVar, "must have more than 1 level"))
    if (all(table(ind) == 1)) 
        stop(paste(groupVar, "must have more than 1 observation by level"))
    p <- ncol(x)
    q1 <- ncol(z)
    if ((sum(is.na(x)) + sum(is.na(z)) + sum(is.na(y)) + sum(is.na(ind))) > 
        0) 
        stop("NAs not allowed")
    if (!is.null(timeVar) && sum(is.na(data[, timeVar]))) 
        stop("NAs not allowed")
    if (distr == "ssl") 
        distr <- "ss"
    if (!(distr %in% c("sn", "st", "ss", "scn"))) 
        stop("Accepted distributions: sn, st, ssl, scn")
    if ((!is.null(control$lb)) && distr != "sn") 
        if ((distr == "st" && (control$lb <= 2)) || (distr =="ss" && (control$lb <= 1)))  stop("Invalid lb")
    if (is.null(control$lb) && distr != "sn") control$lb = ifelse(distr == "scn", rep(0.01, 2), ifelse(distr =="st", 2.01, 1.01))
    if (is.null(control$lu) && distr != "sn") control$lu = ifelse(distr == "scn", rep(0.99, 2), ifelse(distr =="st", 100, 50))
    if (depStruct == "ARp" && !is.null(timeVar) && ((sum(!is.wholenumber(data[, 
        timeVar])) > 0) || (sum(data[, timeVar] <= 0) > 0))) 
        stop("timeVar must contain positive integer numbers when using ARp dependency")
    if (depStruct == "ARp" && !is.null(timeVar)) 
        if (min(data[, timeVar]) != 1) 
            warning("consider using a transformation such that timeVar starts at 1")
    if (depStruct == "CI") 
        depStruct = "UNC"
    if (!(depStruct %in% c("UNC", "ARp", "CS", "DEC", "CAR1"))) 
        stop("accepted depStruct: UNC, ARp, CS, DEC or CAR1")
    if (!(covRandom %in% c("pdSymm", "pdDiag"))) 
        stop("accepted covRandom: pdSymm or pdDiag")
    diagD <- covRandom == "pdDiag"
    if (q1 == 1) 
        diagD = FALSE
    if (missing(skewind)) 
        skewind <- rep(1, q1)
    if (!all(skewind %in% c(0, 1))) 
        stop("skewind must be a vector containing 0s and/or 1s of length equal to the number of random effects")
    if (length(skewind) != q1) 
        stop("skewind must be a vector containing 0s and/or 1s of length equal to the number of random effects")
    if (all(skewind == 0)) 
        stop("for the symmetrical model please use the function smn.lmm")
    if (is.null(control$parallelphi)) 
        control$parallelphi <- ifelse(m > 30, TRUE, FALSE)
    if (depStruct == "UNC") 
        control$parallelphi <- FALSE
    if (is.null(control$parallelnu)) {
        if (distr == "st" || distr == "sn") 
            control$parallelnu <- FALSE
        else if (distr == "ss") 
            control$parallelnu <- ifelse(m > 30, TRUE, FALSE)
        else control$parallelnu <- ifelse(m > 50, TRUE, FALSE)
    }
    if (distr == "sn") 
        control$parallelnu <- FALSE
    if (is.null(control$ncores)) 
        if (control$parallelnu || control$parallelphi) 
            control$ncores <- max(parallel::detectCores() - 1, 
                1, na.rm = TRUE)
    if (is.null(control$initialValues$beta) || is.null(control$initialValues$sigma2) || 
        is.null(control$initialValues$lambda) || is.null(control$initialValues$D)) {
        lmefit = try(lme(formFixed, random = formula(paste("~", 
            as.character(formRandom)[length(formRandom)], "|", 
            "ind")), data = data), silent = T)
        if (is(lmefit, "try-error")) {
            lmefit = try(lme(formFixed, random = ~1 | ind, data = data), 
                silent = TRUE)
            if (is(lmefit, "try-error")) {
                stop("error in calculating initial values")
            }
            else {
                lambdainit <- rep(1, q1) * as.numeric(skewness(random.effects(lmefit)))
                D1init <- diag(q1) * as.numeric(var(random.effects(lmefit)))
            }
        }
        else {
            lambdainit <- skewness(unlist(random.effects(lmefit))) 
            D1init <- (var(random.effects(lmefit)))
        }
    }
    if (!is.null(control$initialValues$beta)) {
        beta1 <- control$initialValues$beta
    }
    else beta1 <- as.numeric(lmefit$coefficients$fixed)
    if (!is.null(control$initialValues$sigma2)) {
        sigmae <- control$initialValues$sigma2
    }
    else sigmae <- as.numeric(lmefit$sigma^2)
    if (!is.null(control$initialValues$D)) {
        D1 <- control$initialValues$D
    }
    else D1 <- D1init
    if (!is.null(control$initialValues$lambda)) {
        lambda <- control$initialValues$lambda
    }
    else lambda <- lambdainit
    lambda <- lambda*skewind
    ####################################### Clecio #################################
    alphas=control$alphas
    nknots=control$nknots
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }

gammas<-list()
for (j in 1:qnl){
     gammas[[j]]=solve(t(Nm[[j]])%*%Nm[[j]]+alphas[j]*Km[[j]])%*%t(Nm[[j]])%*%(y-x%*%beta1)
     
     }
gamma1=unlist(gammas) 
    
################################################################################
    if (length(D1) == 1 && !is.matrix(D1)) 
        D1 = as.matrix(D1)
    if (length(beta1) != p) 
        stop("wrong dimension of beta")
    if (length(lambda) != q1) 
        stop("wrong dimension of lambda")
    if (!is.matrix(D1)) 
        stop("D must be a matrix")
    if ((ncol(D1) != q1) || (nrow(D1) != q1)) 
        stop("wrong dimension of D")
    if (length(sigmae) != 1) 
        stop("wrong dimension of sigma2")
    if (sigmae <= 0) 
        stop("sigma2 must be positive")
    if (depStruct == "ARp") 
        phiAR <- control$initialValues$phi
    if (depStruct == "CS") 
        phiCS <- control$initialValues$phi
    if (depStruct == "DEC") 
        parDEC <- control$initialValues$phi
    if (depStruct == "CAR1") 
        phiCAR1 <- control$initialValues$phi
    nu = control$initialValues$nu
    if (distr == "st" && is.null(nu)) 
        nu = 10
    if (distr == "ss" && is.null(nu)) 
        nu = 5
    if (distr == "scn" && is.null(nu)) 
        nu = c(0.05, 0.8)
    if (distr == "st" && length(nu) != 1) 
        stop("wrong dimension of nu")
    if (distr == "ss" && length(nu) != 1) 
        stop("wrong dimension of nu")
    if (distr == "scn" && length(nu) != 2) 
        stop("wrong dimension of nu")

    if (depStruct == "UNC") 
        obj.out <- DAAREM.SkewUNC(formFixed,formFixedNL, formRandom, data, 
            groupVar, distr, beta1, sigmae, D1, lambda, nu, gammas, lb = control$lb, 
            lu = control$lu, skewind = skewind, diagD = diagD, 
            precisao = control$tol, informa = control$calc.se, 
            max.iter = control$max.iter, showiter = !control$quiet, 
            showerroriter = (!control$quiet) && control$showCriterium, 
            algorithm = control$algorithm, control.daarem = control$control.daarem, 
            parallelnu = control$parallelnu, ncores = control$ncores,alphas,nknots)
    if (depStruct == "ARp") 
            obj.out <- DAAREM.SkewAR(formFixed, formFixedNL, formRandom, data, 
            groupVar, pAR, timeVar, distr, beta1, sigmae, phiAR, 
            D1, lambda, nu, gammas, lb = control$lb, lu = control$lu, 
            skewind = skewind, diagD = diagD, precisao = control$tol, 
            informa = control$calc.se, max.iter = control$max.iter, 
            showiter = !control$quiet, showerroriter = (!control$quiet) && 
                control$showCriterium, algorithm = control$algorithm, 
            control.daarem = control$control.daarem, parallelphi = control$parallelphi, 
            parallelnu = control$parallelnu, ncores = control$ncores,alphas,nknots)
    if (depStruct == "CS") 
        obj.out <- DAAREM.SkewCS(formFixed, formFixedNL, formRandom, data, 
            groupVar, distr, beta1, sigmae, phiCS, D1, lambda, 
            nu, gammas, lb = control$lb, lu = control$lu, skewind = skewind, 
            diagD = diagD, precisao = control$tol, informa = control$calc.se, 
            max.iter = control$max.iter, showiter = !control$quiet, 
            showerroriter = (!control$quiet) && control$showCriterium, 
            algorithm = control$algorithm, control.daarem = control$control.daarem, 
            parallelphi = control$parallelphi, parallelnu = control$parallelnu, 
            ncores = control$ncores,alphas,nknots)
    if (depStruct == "DEC") 
        obj.out <- DAAREM.SkewDEC(formFixed, formFixedNL, formRandom, data, 
            groupVar, timeVar, distr, beta1, sigmae, parDEC, 
            D1, lambda, nu, gammas, lb = control$lb, lu = control$lu, 
            luDEC = control$luDEC, skewind = skewind, diagD = diagD, 
            precisao = control$tol, informa = control$calc.se, 
            max.iter = control$max.iter, showiter = !control$quiet, 
            showerroriter = (!control$quiet) && control$showCriterium, 
            algorithm = control$algorithm, control.daarem = control$control.daarem, 
            parallelphi = control$parallelphi, parallelnu = control$parallelnu, 
            ncores = control$ncores,alphas,nknots)
    if (depStruct == "CAR1") 
        obj.out <- DAAREM.SkewCAR1(formFixed, formFixedNL, formRandom, data, 
            groupVar, timeVar, distr, beta1, sigmae, phiCAR1, 
            D1, lambda, nu, gammas, lb = control$lb, lu = control$lu, 
            skewind = skewind, diagD = diagD, precisao = control$tol, 
            informa = control$calc.se, max.iter = control$max.iter, 
            showiter = !control$quiet, showerroriter = (!control$quiet) && 
                control$showCriterium, algorithm = control$algorithm, 
            control.daarem = control$control.daarem, parallelphi = control$parallelphi, 
            parallelnu = control$parallelnu, ncores = control$ncores,alphas,nknots)
    obj.out$call <- match.call()
    npar <- length(obj.out$theta)
    N <- nrow(data)
    
#    obj.out$criteria$AIC <- 2 * npar - 2 * obj.out$loglik
#    obj.out$criteria$BIC <- log(N) * npar - 2 * obj.out$loglik
    obj.out$criteria$AIC <- 2 * (npar+length(unlist(obj.out$estimates$gammas))) - 2 * obj.out$loglik
    obj.out$criteria$BIC <- log(N) * (npar+length(unlist(obj.out$estimates$gammas))) - 2 * obj.out$loglik
    obj.out$data <- data
    obj.out$formula$formFixed <- formFixed
    obj.out$formula$formRandom <- formRandom
    obj.out$formula$formFixedNL <- formFixedNL
    obj.out$depStruct <- depStruct
    obj.out$covRandom <- covRandom
    if (distr == "ss") 
        distr <- "ssl"
    obj.out$distr <- distr
    obj.out$N <- N
    obj.out$n <- m
    obj.out$groupVar <- groupVar
    obj.out$timeVar <- timeVar
    obj.out$control <- control
    obj.out$diagD <- diagD
    obj.out$skewind <- skewind
    ####################################### Clecio
    # Para calcular Var(Yi) estimada
    y <-data[,all.vars(formFixed)[1]]
    ind <-data[,groupVar]
    data$ind <- ind
    if (is.null(timeVar)) {
       time <- numeric(length = length(ind))
       for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    } else time <- data[,timeVar]
    theta=obj.out$theta
    estimates=obj.out$estimates
    D1=estimates$D
    sigmae2=as.numeric(estimates$sigma2)
    lambda=as.matrix(estimates$lambda)
    delta=lambda/sqrt(1+as.numeric(t(lambda)%*%lambda))
    Delta=matrix.sqrt(D1)%*%delta
    fitted <- numeric(N)
    res_std <- numeric(N)
    ind_levels <- levels(ind)
    gammas=obj.out$estimates$gammas

        c. = -sqrt(2/pi)
        k2=1
        if (distr=="st") {
            nu=as.numeric(theta[length(theta)])
            c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
            k2=nu/2*gamma((nu-2)/2)/gamma(nu/2)}
        if (distr=="ssl") {
            nu=as.numeric(theta[length(theta)])
            c.=-sqrt(2/pi)*nu/(nu-.5)
            k2=nu/(nu-1)}
        if (distr=="scn") {
            l1=length(theta)
            nu=as.vector(theta[(l1-1):l1])
            c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
            k2=1+nu[1]*(nu[2]^(-1)-1)
            }   
        
    aux=1:N
    for (i in seq_along(ind_levels)) {
        auxi <- ind == ind_levels[i]
        seqi <- aux[auxi]
        if (depStruct=="UNC") Ri <-diag(length(seqi))
        if (depStruct=="ARp") Ri <-CovARp(estimates$phi,time[seqi])
        if (depStruct=="CS")  Ri <- CovCS(estimates$phi,length(time[seqi]))
        if (depStruct=="DEC") Ri <- CovDEC(estimates$phi[1],estimates$phi[2],time[seqi])
        if (depStruct=="CAR1") Ri <- CovDEC(estimates$phi,1,time[seqi])         
        Sigmai = sigmae2*Ri
        medg=0
        for (j in 1:qnl){
             Ni=Nm[[j]][seqi,]
             mui=Ni%*%gammas[[j]]
             medg=medg+mui
       }
       xfiti <- matrix(x[seqi, ], ncol = p)
       zfiti <- matrix(z[seqi, ], ncol = q1)
       yi <- as.matrix(y[seqi])
       Gammai=k2*(Sigmai+zfiti%*%D1%*%t(zfiti))-c.^2*zfiti%*%Delta%*%t(Delta)%*%t(zfiti) # Var(Yi)
       res_std[seqi]=solve(matrix.sqrt(Gammai))%*%(yi-xfiti%*%estimates$beta - medg)
       fitted[seqi] <- xfiti %*% obj.out$estimates$beta + medg + zfiti%*%obj.out$random.effects[i,]
    }
    obj.out$res.std <- res_std
    obj.out$fitted <- fitted
    obj.out$Nmatrix <- Nm
    obj.out$Kmatrix <- Km
    obj.out$alphas <- alphas
    obj.out$nknots <- nknots
    class(obj.out) <- c("SMSN", "list")
    obj.out
}



################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################


#SMSN
DAAREM.SkewAR<- function(formFixed,formFixedNL,formRandom,data,groupVar,pAR,timeVar,
                        distr,beta1,sigmae,phiAR,D1,lambda,nu, gammas,lb,lu,
                        skewind,diagD,
                        precisao,informa,max.iter,showiter,showerroriter,
                        algorithm="EM", #"EM" or "DAAREM"
                        parallelphi,parallelnu,ncores,#=detectCores()
                        control.daarem=list(),alphas,nknots){

  ti <- Sys.time()
  data<-as.data.frame(data)
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]

  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
  ng=length(gammas)
  medg=c()
  for (i in 1:ng){
       Ni=Nm[[i]]
       mui=Ni%*%gammas[[i]]
       medg=medg+mui
       }
 
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if ((!is.null(phiAR)) && pAR!=length(phiAR)) stop("initial value from phi must be in agreement with pAR")
  if ((pAR%%1)!=0||pAR==0) stop("pAR must be integer greater than 1")

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)

  if (is.null(phiAR)) {
    lmeAR <- try(lme(formFixed,random=~1|ind,data=data,correlation=corARMA(p=pAR,q=0)),silent=T)# VI de AR soh usando X*Beta
    if (is(lmeAR,"try-error")) piAR =as.numeric(pacf(y-x%*%beta1-medg,lag.max=pAR,plot=F)$acf)
    else {
      phiAR <- capture.output(lmeAR$modelStruct$corStruct)[3]
      phiAR <- as.numeric(strsplit(phiAR, " ")[[1]])
      phiAR <- phiAR[!is.na(phiAR)]
      piAR <- tphitopi(phiAR)
    }
  } else piAR <- tphitopi(phiAR)
  if (any(piAR< -1 | piAR>1)) stop("invalid initial value from phi")

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,piAR,nu)
  ##
  llji <- logveroARpi(y, x, z,Nm,Km, time,ind, beta1, sigmae,piAR, D1, lambda, distr, nu, gammas,alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,length(phiAR)*parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovARp","estphit","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovARp","estphit","logveroARpi","logveroAR","matrix.sqrt",
                          "dmvnorm","ljtAR","ljsAR","ljcnAR"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovARp","estphit","logveroARpi",
                          "logveroAR","matrix.sqrt","dmvnorm","ljtAR","ljsAR","ljcnAR"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewAR,objfn = objfn.skewAR,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewAR,objfn = objfn.skewAR,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  piAR<-EMout$par[(p+2+q2+q1):(p+1+q2+q1+pAR)]
  thetav=EMout$par
  nvnl=length(Km) # numero de var. nao-lineares
  if (distr=="sn") {
      nu <- NULL
      gammav<-thetav[-(1:(p+1+q2+q1+pAR))]
  } 
  if (distr=="st"){
      nu<-thetav[p+1+q2+q1+pAR+1]
      gammav<-thetav[-(1:(p+1+q2+q1+pAR+1))]
      }
  if (distr=="ss"){
      nu<-thetav[p+1+q2+q1+pAR+1]
      gammav<-thetav[-(1:(p+1+q2+q1+pAR+1))]
      }
  if (distr=="scn"){
      nu<-thetav[(p+1+q2+q1+pAR+1):(p+1+q2+q1+pAR+2)]
      gammav<-thetav[-(1:(p+1+q2+q1+pAR+2))]
      }  
#  gammas=list()
#  for (j in 1:qnl){                        # Fiz de outro jeito
#       nj=ncol(Km[[j]])
#       gammas[[j]]=gammav[indg==j]
#       }

auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }  

  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjAR,y=y, x=x, z=z,Nm=Nm, time=time, beta1=beta1, Gammab=Gammab, Deltab=Deltab, sigmae=sigmae,piAR=piAR,zeta=zeta, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y,time=time, x=x, z=z,Nm=Nm, beta1=beta1, Gammab=Gammab,
               Deltab=Deltab, sigmae=sigmae,phi=estphit(piAR),depStruct="ARp", distr=distr,nu=nu,gammas=gammas,simplify = TRUE)
  phiAR<-estphit(piAR)
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiAR,dd,lambda[skewind==1],nu)
  #theta <- c(beta1,sigmae,phiAR,dd,lambda,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),
                                   names_dd, paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),
                        names_dd, paste0("lambda",(1:q1)[skewind==1]),
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=phiAR,dsqrt=dd,D=D1,lambda=as.numeric(lambda),gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  MI=InfmatrixAR(y,x,z,Nm,Km,time,ind,beta1,sigmae,phiAR,D1,lambda,distr = distr,
                           nu = nu,diagD=diagD,skewind=skewind,gammas,alphas)
  obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)
    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}


DAAREM.SkewUNC<- function(formFixed,formFixedNL,formRandom,data,groupVar,
                         distr,beta1,sigmae,D1,lambda,nu, gammas,lb,lu,
                         skewind,diagD,
                         precisao,informa,max.iter,showiter,showerroriter,
                         algorithm="EM", #"EM" or "DAAREM"
                         parallelnu,ncores,
                         control.daarem=list(),alphas,nknots){

ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
# indg=c()
# for (j in 1:length(Km)){
#      aux=rep(j,ncol(Km[[j]]))
#      indg=c(indg,aux)
#      }
          
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #
  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,nu)
  ##
  llji <- logvero(y, x, z,Nm,Km, ind, beta1, sigmae, D1, lambda, distr, nu,gammas,alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu) {
    ncores <- min(ncores,1+2*length(nu))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    clusterExport(cl, c("logvero","matrix.sqrt","dmvnorm","ljt","ljs","ljcn",
                        "n_distinct"),
                    envir=environment())
  }
    if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewUNC,objfn = objfn.skewUNC,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelnu=parallelnu,
                    diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewUNC,objfn = objfn.skewUNC,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  }
 
  if (parallelnu) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]  
  if (distr=="sn") {
    nu <- NULL
  }# else nu<-EMout$par[-(1:(p+1+q2+q1))]
  thetav=EMout$par
  gammav=thetav[-(1:(p+1+q2+q1))]
    if (distr=="st"){
      nu=thetav[p+1+q2+q1+1]
      gammav=thetav[-(1:(p+1+q2+q1+1))]
      }
  if (distr=="ss"){
      nu=thetav[p+1+q2+q1+1]
      gammav=thetav[-(1:(p+1+q2+q1+1))]
      }
  if (distr=="scn"){
      nu=thetav[(p+1+q2+q1+1):(p+1+q2+q1+2)]
      gammav=thetav[-(1:(p+1+q2+q1+2))]
      }

auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }  
  
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emj,y=y, x=x, z=z,Nm=Nm, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae, zeta=zeta, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y, x=x, z=z,Nm=Nm, beta1=beta1, Gammab=Gammab,Deltab=Deltab,
               sigmae=sigmae,depStruct="UNC", distr=distr,nu=nu,gammas=gammas,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control") # class(dd)[1]=='try-error'
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,dd,lambda[skewind==1],nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",names_dd,
                                   paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2",names_dd,
                        paste0("lambda",(1:q1)[skewind==1]),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 dsqrt=dd,D=D1,lambda=as.numeric(lambda),gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  MI=Infmatrix(y,x,z,Nm,Km,ind,beta1,sigmae,D1,lambda,distr = distr,
                           nu = nu,diagD=diagD,skewind=skewind,gammas,alphas)
  obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)

    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}


DAAREM.SkewCS<- function(formFixed,formFixedNL,formRandom,data,groupVar,
                         distr,beta1,sigmae,phiCS,D1,lambda,nu, gammas,lb,lu,
                         skewind,diagD,
                         precisao,informa,max.iter,showiter,showerroriter,
                         algorithm="EM", #"EM" or "DAAREM"
                         parallelphi,parallelnu,ncores,
                         control.daarem=list(),alphas,nknots){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }

  ng=length(gammas)
  medg=0
  for (i in 1:ng){
       Ni=Nm[[i]]
       mui=Ni%*%gammas[[i]]
       medg=medg+mui
       }
 
          
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #
  if (!is.null(phiCS) && length(phiCS)!=1) stop ("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCS)) if (phiCS<=0 | phiCS>=1) stop("0<initialValue$phi<1 needed")
  #
  if (is.null(phiCS)) {
    phiCS <- abs(as.numeric(pacf(y-x%*%beta1-medg,lag.max=1,plot=F)$acf))
  }
  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiCS,nu)
  ##
  llji <- logveroCS(y, x, z,Nm,Km,ind, beta1, sigmae,phiCS, D1, lambda, distr, nu,gammas,alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovCS","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovCS","logveroCS","matrix.sqrt",
                          "dmvnorm","ljtCS","ljsCS","ljcnCS"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovCS","logveroCS","matrix.sqrt",
                          "dmvnorm","ljtCS","ljsCS","ljcnCS"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewCS,objfn = objfn.skewCS,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewCS,objfn = objfn.skewCS,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,
                    diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  phiCS<-EMout$par[(p+2+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  }# else nu<-EMout$par[-(1:(p+2+q2+q1))]
  thetav=EMout$par
  gammav=thetav[-(1:(p+2+q2+q1))]
  if (distr=="st"){
      nu=thetav[p+2+q2+q1+1]
      gammav=thetav[-(1:(p+2+q2+q1+1))]
      }
  if (distr=="ss"){
      nu=thetav[p+2+q2+q1+1]
      gammav=thetav[-(1:(p+2+q2+q1+1))]
      }
  if (distr=="scn"){
      nu=thetav[(p+2+q2+q1+1):(p+2+q2+q1+2)]
      gammav=thetav[-(1:(p+2+q2+q1+2))]
      }

auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }  
  
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjCS,y=y, x=x, z=z,Nm=Nm, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,phiCS=phiCS,zeta=zeta, distr=distr,
                             nu=nu,gammas=gammas,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y, x=x, z=z,Nm=Nm, beta1=beta1, Gammab=Gammab,Deltab=Deltab,
               sigmae=sigmae,depStruct="CS", phi=phiCS,distr=distr,nu=nu,gammas=gammas,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control") #class(dd)[1]=='try-error'
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiCS,dd,lambda[skewind==1],nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCS",
                                   names_dd, paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2","phiCS",names_dd,
                        paste0("lambda",(1:q1)[skewind==1]),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae, phi=phiCS,
                                 dsqrt=dd,D=D1,lambda=as.numeric(lambda),gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  MI=InfmatrixCS(y,x,z,Nm,Km,ind,beta1,sigmae,phiCS,D1,lambda,distr = distr,
                           nu = nu,diagD=diagD,skewind=skewind,gammas,alphas)
   obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)

    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}


DAAREM.SkewDEC<- function(formFixed,formFixedNL,formRandom,data,groupVar,timeVar,
                         distr,beta1,sigmae,parDEC,D1,lambda,nu, gammas,lb,lu,luDEC,
                         skewind,diagD,
                         precisao,informa,max.iter,showiter,showerroriter,
                         algorithm="EM", #"EM" or "DAAREM"
                         parallelphi,parallelnu,ncores,
                         control.daarem=list(),alphas,nknots){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
 
#  indg=c()
# for (j in 1:length(Km)){
#      aux=rep(j,ncol(Km[[j]]))
#      indg=c(indg,aux)
#      }
      
         
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(parDEC)) {
    if (length(parDEC)!=2) stop ("initial value from phi should have length 2 or NULL")
    if (parDEC[1]<=0||parDEC[1]>=1) stop("invalid initial value from phi1")
    if (parDEC[2]<=0) stop("invalid initial value from phi2")
    if (parDEC[2]>= luDEC) stop("initial value from phi2 must be smaller than luDEC")
  }

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)

  if (is.null(parDEC)) {
    #cat("calculating initial values for DEC... \n")
    thetat<- seq(0.1,luDEC,by=.1)
    phit <- seq(0.1,.9,by=.05)
    vect <-merge(phit,thetat,all=T)
    logveroDECv<-function(phitheta){logveroDEC(y, x, z,Nm,Km,time,ind, beta1=beta1, sigmae=sigmae,
    phiDEC=phitheta[1],thetaDEC=phitheta[2], D1=D1,lambda=lambda, distr=distr, nu=nu,gammas,alphas)}
    logverovec <- apply(vect,1,logveroDECv)
    parDEC <- as.numeric(vect[which.max(logverovec),])
  }
  #phiDEC=parDEC[1]
  #thetaDEC=parDEC[2]
  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,parDEC,nu)
  ##
  llji <- logveroDEC(y, x, z,Nm,Km, time,ind, beta1, sigmae,parDEC[1],parDEC[2], D1, lambda, distr, nu,gammas,alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,2*parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","logveroDEC","matrix.sqrt",
                          "dmvnorm","ljtDEC","ljsDEC","ljcnDEC"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovDEC","logveroDEC","matrix.sqrt",
                          "dmvnorm","ljtDEC","ljsDEC","ljcnDEC"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewDEC,objfn = objfn.skewDEC,
                    y=y,x=x,z=z,Nm,Km,time=time,ind=ind,distr=distr,lb=lb,lu=lu,luDEC=luDEC,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewDEC,objfn = objfn.skewDEC,
                    y=y,x=x,z=z,Nm=Nm,time=time,ind=ind,distr=distr,lb=lb,lu=lu,luDEC=luDEC,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,
                    diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  phiDEC<-EMout$par[(p+2+q2+q1)]
  thetaDEC<-EMout$par[(p+3+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  }# else nu<-EMout$par[-(1:(p+3+q2+q1))]
  thetav=EMout$par
  gammav=thetav[-(1:(p+3+q2+q1))]
  if (distr=="st"){
      nu=thetav[p+3+q2+q1+1]
      gammav=thetav[-(1:(p+3+q2+q1+1))]
      }
  if (distr=="ss"){
      nu=thetav[p+3+q2+q1+1]
      gammav=thetav[-(1:(p+3+q2+q1+1))]
      }
  if (distr=="scn"){
      nu=thetav[(p+3+q2+q1+1):(p+3+q2+q1+2)]
      gammav=thetav[-(1:(p+3+q2+q1+2))]
      }  

auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }  
  
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDEC,y=y, x=x, z=z,Nm=Nm, time=time, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC,zeta=zeta, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y,time=time, x=x, z=z,Nm=Nm, beta1=beta1, Gammab=Gammab,
               Deltab=Deltab, sigmae=sigmae,phi=c(phiDEC,thetaDEC),depStruct="DEC", distr=distr,nu=nu,gammas=gammas,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control") #class(dd)[1]=='try-error'
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiDEC,thetaDEC,dd,lambda[skewind==1],nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phi1DEC","phi2DEC",
                                   names_dd, paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2","phi1DEC","phi2DEC",
                        names_dd, paste0("lambda",(1:q1)[skewind==1]),
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=c(phiDEC,thetaDEC),dsqrt=dd,D=D1,lambda=as.numeric(lambda),gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  MI=InfmatrixDEC(y,x,z,Nm,Km,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,lambda,distr = distr,
                           nu = nu,diagD=diagD,skewind=skewind,gammas,alphas)
   obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)

    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}


DAAREM.SkewCAR1 <- function(formFixed,formFixedNL,formRandom,data,groupVar,timeVar,
                          distr,beta1,sigmae,phiCAR1,D1,lambda,nu, gammas,lb,lu,
                          skewind,diagD, precisao,informa,max.iter,showiter,showerroriter,
                          algorithm="EM", #"EM" or "DAAREM"
                          parallelphi,parallelnu,ncores,
                          control.daarem=list(),alphas,nknots){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }

# indg=c()
# for (j in 1:length(Km)){
#      aux=rep(j,ncol(Km[[j]]))
#      indg=c(indg,aux)
#      }
  ng=length(gammas)
  medg=0
  for (i in 1:ng){
       Ni=Nm[[i]]
       mui=Ni%*%gammas[[i]]
       medg=medg+mui
       }
          
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(phiCAR1) && length(phiCAR1)!=1) stop("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCAR1)) if (phiCAR1>=1 || phiCAR1<=0) stop ("0<initialValue$phi<1 needed")

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)

  if (is.null(phiCAR1)) {
    lmeCAR = try(lme(formFixed,random=~1|ind,data=data,correlation=corCAR1(form = ~time)),silent=T)
    if (is(lmeCAR,"try-error")) phiDEC =abs(as.numeric(pacf(y-x%*%beta1-medg,lag.max=1,plot=F)$acf))
    else {
      phiDEC = capture.output(lmeCAR$modelStruct$corStruct)[3]
      phiDEC = as.numeric(strsplit(phiDEC, " ")[[1]])
    }
  } else phiDEC <- phiCAR1

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiDEC,nu)
  ##
  llji <- logveroCAR1(y, x, z,Nm,Km, time,ind, beta1, sigmae,phiDEC, D1, lambda, distr, nu, gammas, alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","logveroCAR1","matrix.sqrt",
                          "dmvnorm","ljtCAR1","ljsCAR1","ljcnCAR1"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovDEC","logveroCAR1","matrix.sqrt",
                          "dmvnorm","ljtCAR1","ljsCAR1","ljcnCAR1"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewCAR1,objfn = objfn.skewCAR1,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewCAR1,objfn = objfn.skewCAR1,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,
                    diagD=diagD,skewind=skewind,gammas=gammas,alphas=alphas)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")

  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  phiDEC<-EMout$par[(p+2+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  }# else nu<-EMout$par[-(1:(p+2+q2+q1))]
  thetav=EMout$par
  gammav=thetav[-(1:(p+2+q2+q1))]
  if (distr=="st"){
      nu=thetav[p+2+q2+q1+1]
      gammav=thetav[-(1:(p+2+q2+q1+1))]
      }
  if (distr=="ss"){
      nu=thetav[p+2+q2+q1+1]
      gammav=thetav[-(1:(p+2+q2+q1+1))]
      }
  if (distr=="scn"){
      nu=thetav[(p+2+q2+q1+1):(p+2+q2+q1+2)]
      gammav=thetav[-(1:(p+2+q2+q1+2))]
      }  

auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }  
  
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDEC,y=y, x=x, z=z,Nm=Nm, time=time, beta1=beta1, Gammab=Gammab, Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=1,zeta=zeta, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y,time=time, x=x, z=z,Nm=Nm, beta1=beta1, Gammab=Gammab,
               Deltab=Deltab, sigmae=sigmae,phi=c(phiDEC,1),depStruct="DEC", distr=distr,nu=nu,gammas=gammas,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control") #class(dd)[1]=='try-error'
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiDEC,dd,lambda[skewind==1],nu)

  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCAR1",
                                   names_dd, paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2","phiCAR1",
                        names_dd, paste0("lambda",(1:q1)[skewind==1]),
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=phiDEC,dsqrt=dd,D=D1,lambda=as.numeric(lambda),gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  MI=InfmatrixCAR1(y,x,z,Nm,Km,time,ind,beta1,sigmae,phiDEC,D1,lambda,distr = distr,
                           nu = nu,diagD=diagD,skewind=skewind,gammas,alphas)
   obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)

    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
###################SMN



smn.lmm  <- function(data, formFixed, formFixedNL, groupVar, formRandom = ~1, depStruct = "UNC", 
timeVar = NULL, distr = "n", covRandom = "pdSymm", pAR = 1, control = lmmControl())    # data tem de ser formato dataframe
{

  if (!is(formFixed,"formula")) stop("formFixed must be a formula")
  if (!is(formRandom,"formula")) stop("formRandom must be a formula")
  if (!inherits(control,"lmmControl")) stop("control must be a list generated with lmmControl()")
  #
  if (!is.character(groupVar)) stop("groupVar must be a character containing the name of the grouping variable in data")
  if (!is.null(timeVar)&&!is.character(timeVar)) stop("timeVar must be a character containing the name of the time variable in data")
  if (length(formFixed)!=3) stop("formFixed must be a two-sided linear formula object")
  if (!is.data.frame(data)) stop("data must be a data.frame")
  #if (is_tibble(data)) data=as.data.frame(data)
  if (length(class(data))>1) data=as.data.frame(data)
  vars_used<-unique(c(all.vars(formFixed),all.vars(formRandom),groupVar,timeVar))
  vars_miss <- which(!(vars_used %in% names(data)))
  if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in data"))
  data = data[,vars_used]
  #
  #data <- data[order(data[,groupVar]),]
  if (!is.factor(data[,groupVar])) data[,groupVar]<-haven::as_factor(data[,groupVar])
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <-data[,groupVar]
  m<-nlevels(ind)
  if (m<=1) stop(paste(groupVar,"must have more than 1 level"))
  if (all(table(ind)==1)) stop(paste(groupVar,"must have more than 1 observation by level"))
  p<-ncol(x)
  q1<-ncol(z)
  if ((sum(is.na(x))+sum(is.na(z))+sum(is.na(y))+sum(is.na(ind)))>0) stop ("NAs not allowed")
  if (!is.null(timeVar) && sum(is.na(data[,timeVar]))) stop ("NAs not allowed")
  #
  if (!(distr %in% c("norm","t","sl","cn"))) stop("Accepted distributions: norm, t, sl, cn")
  if ((!is.null(control$lb))&&distr!="norm") if((distr=="t"&&(control$lb<=2))||(distr=="sl"&&(control$lb<=1))) stop("Invalid lb")
  if (is.null(control$lb)&&distr!="norm") control$lb = ifelse(distr=="cn",rep(.01,2),ifelse(distr=="t",2.01,1.01))
  if (is.null(control$lu)&&distr!="norm") control$lu = ifelse(distr=="cn",rep(.99,2),ifelse(distr=="t",100,50))
  #
  if (depStruct=="ARp" && !is.null(timeVar) &&
      ((sum(!is.wholenumber(data[,timeVar]))>0)||(sum(data[,timeVar]<=0)>0))) stop("timeVar must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp" && !is.null(timeVar)) if (min(data[,timeVar])!=1) warning("consider using a transformation such that timeVar starts at 1")
  if (depStruct=="CI") depStruct = "UNC"
  if (!(depStruct %in% c("UNC","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or CAR1")
  #
  if (!(covRandom %in% c('pdSymm','pdDiag'))) stop("accepted covRandom: pdSymm or pdDiag")
  diagD <- covRandom=='pdDiag'
  if (q1==1) diagD=FALSE
  if (is.null(control$parallelphi)) control$parallelphi <- ifelse(m>30,TRUE,FALSE)
  if (depStruct=="UNC") control$parallelphi <- FALSE
  if (is.null(control$parallelnu)) {
    if (distr=='t'||distr=='norm') control$parallelnu <- FALSE
    else if (distr=='sl') control$parallelnu <- ifelse(m>30,TRUE,FALSE)
    else control$parallelnu <- ifelse(m>50,TRUE,FALSE)
  }
  if (distr=='norm') control$parallelnu <- FALSE
  if (is.null(control$ncores)) if (control$parallelnu||control$parallelphi) control$ncores <- max(parallel::detectCores() - 1, 1, na.rm = TRUE)
  
  if (is.null(control$initialValues$beta)||is.null(control$initialValues$sigma2)||
      is.null(control$initialValues$D)) {
    lmefit = try(lme(formFixed,random=formula(paste('~',as.character(formRandom)[length(formRandom)],
                                                    '|',"ind")),data=data),silent=T)

    if (is(lmefit,"try-error")) {
      lmefit = try(lme(formFixed,random=~1|ind,data=data),silent=TRUE)
      if (is(lmefit,"try-error")) {
        stop("error in calculating initial values")
      } else {
        D1init <- diag(q1)*as.numeric(var(random.effects(lmefit)))
      }
    } else {
      D1init <- (var(random.effects(lmefit)))
    }
  }

  if (!is.null(control$initialValues$beta)) {
    beta1 <- control$initialValues$beta
  } else beta1 <- as.numeric(lmefit$coefficients$fixed)
  if (!is.null(control$initialValues$sigma2)) {
    sigmae <- control$initialValues$sigma2
  } else sigmae <- as.numeric(lmefit$sigma^2)
  if (!is.null(control$initialValues$D)) {
    D1 <- control$initialValues$D
  } else D1 <- D1init 

####################################### Clecio #################################
    alphas=control$alphas
    nknots=control$nknots
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)             
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
gammas<-list()
for (j in 1:qnl){
     gammas[[j]]=solve(t(Nm[[j]])%*%Nm[[j]]+alphas[j]*Km[[j]])%*%t(Nm[[j]])%*%(y-x%*%beta1)
     
     }
gamma1=unlist(gammas)  
################################################################################

  #
  if (length(D1)==1 && !is.matrix(D1)) D1=as.matrix(D1)
  #
  if (length(beta1)!=p) stop ("wrong dimension of beta")
  if (!is.matrix(D1)) stop("D must be a matrix")
  if ((ncol(D1)!=q1)|(nrow(D1)!=q1)) stop ("wrong dimension of D")
  if (length(sigmae)!=1) stop ("wrong dimension of sigma2")
  if (sigmae<=0) stop("sigma2 must be positive")
  #
  if (depStruct=="ARp") phiAR<- control$initialValues$phi
  if (depStruct=="CS") phiCS<- control$initialValues$phi
  if (depStruct=="DEC") parDEC<- control$initialValues$phi
  if (depStruct=="CAR1") phiCAR1<- control$initialValues$phi
  #
  nu = control$initialValues$nu
  #
  if (distr=="t"&&is.null(nu)) nu=10
  if (distr=="sl"&&is.null(nu)) nu=5
  if (distr=="cn"&&is.null(nu)) nu=c(.05,.8)
  #
  if (distr=="t"&&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="sl"&&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="cn"&&length(nu)!=2) stop ("wrong dimension of nu")
  #
  if (distr=="norm") distrs="sn"
  if (distr=="t") distrs="st"
  if (distr=="sl") distrs="ss"
  if (distr=="cn") distrs="scn"

  ###
    if (depStruct == "UNC") 
        obj.out <- DAAREM.UNC(formFixed,formFixedNL, formRandom, data, 
            groupVar, distr=distrs, beta1, sigmae, D1, nu, gammas, lb = control$lb, 
            lu = control$lu, diagD = diagD, precisao = control$tol, informa = control$calc.se, 
            max.iter = control$max.iter, showiter = !control$quiet, 
            showerroriter = (!control$quiet) && control$showCriterium, 
            algorithm = control$algorithm, control.daarem = control$control.daarem, 
            parallelnu = control$parallelnu, ncores = control$ncores,alphas=alphas,nknots=nknots)
    if (depStruct == "ARp") 
        obj.out <- DAAREM.AR(formFixed, formFixedNL, formRandom, data, 
            groupVar, pAR, timeVar, distr=distrs, beta1, sigmae, phiAR, 
            D1, nu, gammas, lb = control$lb, lu = control$lu, 
            diagD = diagD, precisao = control$tol, 
            informa = control$calc.se, max.iter = control$max.iter, 
            showiter = !control$quiet, showerroriter = (!control$quiet) && 
                control$showCriterium, algorithm = control$algorithm, 
            control.daarem = control$control.daarem, parallelphi = control$parallelphi, 
            parallelnu = control$parallelnu, ncores = control$ncores,alphas=alphas,nknots=nknots)
    if (depStruct == "CS") 
        obj.out <- DAAREM.CS(formFixed, formFixedNL, formRandom, data, 
            groupVar, distr=distrs, beta1, sigmae, phiCS, D1,nu, gammas, lb = control$lb, lu = control$lu, 
            diagD = diagD, precisao = control$tol, informa = control$calc.se, 
            max.iter = control$max.iter, showiter = !control$quiet, 
            showerroriter = (!control$quiet) && control$showCriterium, 
            algorithm = control$algorithm, control.daarem = control$control.daarem, 
            parallelphi = control$parallelphi, parallelnu = control$parallelnu, 
            ncores = control$ncores,alphas=alphas,nknots=nknots)
    if (depStruct == "DEC") 
        obj.out <- DAAREM.DEC(formFixed, formFixedNL, formRandom, data, 
            groupVar, timeVar, distr=distrs, beta1, sigmae, parDEC, 
            D1, nu, gammas, lb = control$lb, lu = control$lu, 
            luDEC = control$luDEC, diagD = diagD, 
            precisao = control$tol, informa = control$calc.se, 
            max.iter = control$max.iter, showiter = !control$quiet, 
            showerroriter = (!control$quiet) && control$showCriterium, 
            algorithm = control$algorithm, control.daarem = control$control.daarem, 
            parallelphi = control$parallelphi, parallelnu = control$parallelnu, 
            ncores = control$ncores,alphas=alphas,nknots=nknots)
    if (depStruct == "CAR1") 
        obj.out <- DAAREM.CAR1(formFixed, formFixedNL, formRandom, data, 
            groupVar, timeVar, distr=distrs, beta1, sigmae, phiCAR1, 
            D1, nu, gammas, lb = control$lb, lu = control$lu, 
            diagD = diagD, precisao = control$tol, 
            informa = control$calc.se, max.iter = control$max.iter, 
            showiter = !control$quiet, showerroriter = (!control$quiet) && 
                control$showCriterium, algorithm = control$algorithm, 
            control.daarem = control$control.daarem, parallelphi = control$parallelphi, 
            parallelnu = control$parallelnu, ncores = control$ncores,alphas=alphas,nknots=nknots)

    obj.out$call <- match.call()  
    npar <- length(obj.out$theta)
    N <- nrow(data)
    obj.out$criteria$AIC <- 2 * (npar+length(unlist(obj.out$estimates$gammas))) - 2 * obj.out$loglik
    obj.out$criteria$BIC <- log(N) * (npar+length(unlist(obj.out$estimates$gammas))) - 2 * obj.out$loglik
    obj.out$data <- data
    obj.out$formula$formFixed <- formFixed
    obj.out$formula$formRandom <- formRandom
    obj.out$formula$formFixedNl <- formFixedNL
    obj.out$depStruct <- depStruct
    obj.out$covRandom <- covRandom
    obj.out$distr <- distr
    obj.out$N <- N
    obj.out$n <- m
    obj.out$groupVar <- groupVar
    obj.out$timeVar <- timeVar
    obj.out$control <- control
    obj.out$diagD <- diagD
    
    ####################################### Clecio
    # Para calcular Var(Yi) estimada
    y <-data[,all.vars(formFixed)[1]]
    ind <-data[,groupVar]
    data$ind <- ind
    if (is.null(timeVar)) {
       time <- numeric(length = length(ind))
       for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    } else time <- data[,timeVar]
    theta=obj.out$theta
    estimates=obj.out$estimates
    D1=estimates$D
    sigmae2=as.numeric(estimates$sigma2)
    gammas=obj.out$estimates$gammas
    fitted <- numeric(N)
    res_std <- numeric(N)
    ind_levels <- levels(ind)
    
    c. = -sqrt(2/pi)
    k2=1
    if (distr=="t") {
            nu=as.numeric(theta[length(theta)])
            c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
            k2=nu/2*gamma((nu-2)/2)/gamma(nu/2)}
    if (distr=="sl") {
            nu=as.numeric(theta[length(theta)])
            c.=-sqrt(2/pi)*nu/(nu-.5)
            k2=nu/(nu-1)}
    if (distr=="cn") {
            l1=length(theta)
            nu=as.vector(theta[(l1-1):l1])
            c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
            k2=1+nu[1]*(nu[2]^(-1)-1)
    } 

    for (i in seq_along(ind_levels)) {
        seqi <- ind[ind == ind_levels[i]]      
        if (depStruct=="UNC") Ri <-diag(length(seqi))
        if (depStruct=="ARp") Ri <-CovARp(estimates$phi,time[seqi])
        if (depStruct=="CS")  Ri <- CovCS(estimates$phi,length(time[seqi]))
        if (depStruct=="DEC") Ri <- CovDEC(estimates$phi[1],estimates$phi[2],time[seqi])
        if (depStruct=="CAR1") Ri <- CovDEC(estimates$phi,1,time[seqi])         
        Sigmai = sigmae2*Ri
        medg=0
        for (j in 1:qnl){
             Ni=Nm[[j]][seqi,]
             mui=Ni%*%gammas[[j]]
             medg=medg+mui
       }
       xfiti <- matrix(x[seqi, ], ncol = p)
       zfiti <- matrix(z[seqi, ], ncol = q1)
       yi <- as.matrix(y[seqi])
       fitted[seqi] <- xfiti %*% obj.out$estimates$beta + zfiti %*% obj.out$random.effects[i, ] + medg
       Gammai=k2*(Sigmai+zfiti%*%D1%*%t(zfiti)) # Var(Yi)
       #res_std[seqi]=solve(matrix.sqrt(Gammai))%*%(yi-xfiti%*%estimates$beta - medg)  # ev(Gammai)<0 para Normal, AR
    }
    obj.out$fitted <- fitted
    obj.out$res.std <- res_std
    obj.out$Nmatrix <- Nm
    obj.out$Kmatrix <- Km
    obj.out$alphas <- alphas
    obj.out$nknots <- nknots
    class(obj.out) <- c("SMN", "list")
    obj.out
}




DAAREM.AR<- function(formFixed,formFixedNL,formRandom,data,groupVar,pAR,timeVar,
                     distr,beta1,sigmae,phiAR,D1,nu, gammas,lb,lu,diagD,
                     precisao,informa,max.iter,showiter,showerroriter,
                     algorithm="EM", #"EM" or "DAAREM"
                     parallelphi,parallelnu,ncores,
                     control.daarem=list(),alphas,nknots){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
# indg=c()
# for (j in 1:length(Km)){
#      aux=rep(j,ncol(Km[[j]]))
#      indg=c(indg,aux)
#      }
  ng=length(gammas)
  medg=0
  for (i in 1:ng){
       Ni=Nm[[i]]
       mui=Ni%*%gammas[[i]]
       medg=medg+mui
       }
          
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if ((!is.null(phiAR)) && pAR!=length(phiAR)) stop("initial value from phi must be in agreement with pAR")
  if ((pAR%%1)!=0||pAR==0) stop("pAR must be integer greater than 1")

  if (is.null(phiAR)) {
    lmeAR = try(lme(formFixed,random=~1|ind,data=data,correlation=corARMA(p=pAR,q=0)),silent=T)
    if (is(lmeAR,"try-error")) piAR =as.numeric(pacf(y-x%*%beta1-medg,lag.max=pAR,plot=F)$acf)
    else {
      phiAR = capture.output(lmeAR$modelStruct$corStruct)[3]
      phiAR = as.numeric(strsplit(phiAR, " ")[[1]])
      phiAR = phiAR[!is.na(phiAR)]
      piAR = tphitopi(phiAR)
    }
  } else piAR = tphitopi(phiAR)
  if (any(piAR< -1 | piAR>1)) stop("invalid initial value from phi")

  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],piAR,nu)
  ##                  
  llji <- logveroARpis(y, x, z, Nm, Km, time,ind, beta1, sigmae,piAR, D1, distr, nu, gammas,alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")
  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,length(phiAR)*parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovARp","estphit","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovARp","estphit","logveroARpis","logveroARs",
                          "matrix.sqrt","dmvnorm","ljtARs","ljsARs","ljcnARs"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovARp","estphit","logveroARpis",
                          "logveroARs","matrix.sqrt","dmvnorm","ljtARs","ljsARs","ljcnARs"),
                    envir=environment())
    }
  }
  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.AR,objfn = objfn.AR,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.AR,objfn = objfn.AR,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  }
  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  piAR<-EMout$par[(p+2+q2):(p+1+q2+pAR)]
  if (distr=="sn") {
    nu <- NULL
  }# else nu<-EMout$par[-(1:(p+1+q2+pAR))]
  thetav=EMout$par
  gammav<-thetav[-(1:(p+1+q2+pAR))]
  if (distr=="st"){
      nu=thetav[p+1+q2+pAR+1]
      gammav<-thetav[-(1:(p+1+q2+pAR+1))]
      }
  if (distr=="ss"){
      nu=thetav[p+1+q2+pAR+1]
      gammav<-thetav[-(1:(p+1+q2+pAR+1))]
      }
  if (distr=="scn"){
      nu=thetav[(p+1+q2+pAR+1):(p+1+q2+pAR+2)]
      gammav<-thetav[-(1:(p+1+q2+pAR+2))]
      }
 
auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }  

  #
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjARs,y=y, x=x, z=z,Nm=Nm,time=time, beta1=beta1, D1=D1,
                             sigmae=sigmae,piAR=piAR, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y,time=time, x=x, z=z,Nm=Nm, beta1=beta1, D1=D1,
               sigmae=sigmae,phi=estphit(piAR),depStruct="ARp", distr=distr,nu=nu,gammas=gammas,
               simplify = TRUE)
  #cat('\n',D1)
  phiAR<-estphit(piAR)
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiAR,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),
                                   names_dd)
  else names(theta)<- c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),
                        names_dd,paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=phiAR,dsqrt=dd,D=D1,gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  skewind=rep(0,q1)  
                 
  MI=InfmatrixAR(y,x,z,Nm,Km,time,ind,beta1,sigmae,phiAR,D1,lambda=rep(0,q1),distr = distr,
                           nu = nu,diagD=diagD,skewind=rep(0,q1),gammas,alphas)
  obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)
    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}


DAAREM.UNC<- function(formFixed,formFixedNL,formRandom,data,groupVar,
                      distr,beta1,sigmae,D1,nu, gammas,lb,lu,diagD,
                      precisao,informa,max.iter,showiter,showerroriter,
                      algorithm="EM", #"EM" or "DAAREM"
                      parallelnu,ncores,
                      control.daarem=list(),alphas,nknots){

  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]

  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
# indg=c()
# for (j in 1:length(Km)){
#      aux=rep(j,ncol(Km[[j]]))
#      indg=c(indg,aux)
#      }       
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)    # OK
  ##
  llji <- logveros(y, x, z,Nm,Km, ind, beta1, sigmae, D1, distr, nu, gammas,alphas)  # OK

  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")
  if (parallelnu) {
    ncores <- min(ncores,1+2*length(nu))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    clusterExport(cl, c("logveros","matrix.sqrt","dmvnorm","ljts","ljss","ljcns",
                        "n_distinct"),
                  envir=environment())
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.UNC,objfn = objfn.UNC,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.UNC,objfn = objfn.UNC,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,
                    showerroriter = showerroriter,parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  }
  if (parallelnu) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  if (distr=="sn") {
    nu <- NULL
  #} else nu<-EMout$par[-(1:(p+1+q2))]
    } 
  thetav=EMout$par
  gammav<-thetav[-(1:(p+1+q2))]
  if (distr=="st"){
      nu=EMout$par[p+1+q2+1]
      gammav<-thetav[-(1:(p+1+q2+1))]
      }
  if (distr=="ss"){
      nu=EMout$par[p+1+q2+1]
      gammav<-thetav[-(1:(p+1+q2+1))]
      }
  if (distr=="scn"){
      nu=EMout$par[(p+1+q2+1):(p+1+q2+2)]
      gammav<-thetav[-(1:(p+1+q2+2))]
      }

auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }
  #
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)

  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjs,y=y, x=x, z=z,Nm=Nm, beta1=beta1, D1=D1,
                             sigmae=sigmae, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y, x=x, z=z,Nm=Nm, beta1=beta1, D1=D1,
               sigmae=sigmae,depStruct="UNC", distr=distr,nu=nu,gammas=gammas,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",names_dd)
  else names(theta)<- c(colnames(x),"sigma2",names_dd,paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 dsqrt=dd,D=D1,gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi
 
 
 skewind=rep(0,q1)  

  MI=Infmatrix(y,x,z,Nm,Km,ind,beta1,sigmae,D1,lambda=rep(0,q1),distr = distr,
                           nu = nu,diagD=diagD,skewind=rep(0,q1),gammas,alphas)
  obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)
    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)
  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}


DAAREM.CS<- function(formFixed,formFixedNL,formRandom,data,groupVar,
                     distr,beta1,sigmae,phiCS,D1,nu, gammas,lb,lu,diagD,
                     precisao,informa,max.iter,showiter,showerroriter,
                     algorithm="EM", #"EM" or "DAAREM"
                     parallelphi,parallelnu,ncores,
                     control.daarem=list(),alphas,nknots){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
  ng=length(gammas)
  medg=0
  for (i in 1:ng){
       Ni=Nm[[i]]
       mui=Ni%*%gammas[[i]]
       medg=medg+mui
       }         
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(phiCS) && length(phiCS)!=1) stop ("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCS)) if (phiCS<=0 | phiCS>=1) stop("0<initialValue$phi<1 needed")
  #
  if (is.null(phiCS)) {
    phiCS <- abs(as.numeric(pacf(y-x%*%beta1-medg,lag.max=1,plot=F)$acf))
  }
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiCS,nu)
  ##
  llji <- logveroCSs(y, x, z, Nm, Km,ind, beta1, sigmae,phiCS, D1, distr, nu,gammas,alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovCS","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovCS","logveroCSs","matrix.sqrt",
                          "dmvnorm","ljtCSs","ljsCSs","ljcnCSs"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovCS","logveroCSs","matrix.sqrt",
                          "dmvnorm","ljtCSs","ljsCSs","ljcnCSs"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.CS,objfn = objfn.CS,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.CS,objfn = objfn.CS,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  phiCS<-EMout$par[(p+2+q2)]
  if (distr=="sn") {
    nu <- NULL
  }# else nu<-EMout$par[-(1:(p+2+q2))]
  thetav=EMout$par
  gammav<-thetav[-(1:(p+2+q2))]
  if (distr=="st"){
      nu=thetav[p+2+q2+1]
      gammav<-thetav[-(1:(p+2+q2+1))]
      }
  if (distr=="ss"){
      nu=thetav[p+2+q2+1]
      gammav<-thetav[-(1:(p+2+q2+1))]
      }
  if (distr=="scn"){
      nu=thetav[(p+2+q2+1):(p+2+q2+2)]
      gammav<-thetav[-(1:(p+2+q2+2))]
      }
  
auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }  
  
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjCSs,y=y, x=x, z=z,Nm=Nm, beta1=beta1, D1=D1,
                             sigmae=sigmae,phiCS=phiCS, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y, x=x, z=z,Nm=Nm, beta1=beta1, D1=D1,
               sigmae=sigmae,depStruct="CS", phi=phiCS,distr=distr,nu=nu, gammas=gammas,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiCS,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCS",names_dd)
  else names(theta)<- c(colnames(x),"sigma2","phiCS",names_dd,
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
  estimates=list(beta=as.numeric(beta1),sigma2=sigmae, phi=phiCS,
  dsqrt=dd,D=D1,gammas=gammas), uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  MI=InfmatrixCS(y,x,z,Nm,Km,ind,beta1,sigmae,phiCS,D1,lambda=rep(0,q1),distr = distr,
                           nu = nu,diagD=diagD,skewind=rep(0,q1),gammas,alphas)
   obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)
    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}


DAAREM.DEC<- function(formFixed,formFixedNL,formRandom,data,groupVar,timeVar,
                      distr,beta1,sigmae,parDEC,D1,nu, gammas,lb,lu,luDEC,diagD,
                      precisao,informa,max.iter,showiter,showerroriter,
                      algorithm="EM", #"EM" or "DAAREM"
                      parallelphi,parallelnu,ncores,
                      control.daarem=list(),alphas,nknots){
                                        
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
# indg=c()
# for (j in 1:length(Km)){
#      aux=rep(j,ncol(Km[[j]]))
#      indg=c(indg,aux)
#      }      
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(parDEC)) {
    if (length(parDEC)!=2) stop ("initial value from phi should have length 2 or NULL")
    if (parDEC[1]<=0||parDEC[1]>=1) stop("invalid initial value from phi1")
    if (parDEC[2]<=0) stop("invalid initial value from phi2")
    if (parDEC[2]>= luDEC) stop("initial value from phi2 must be smaller than luDEC")
  }

  if (is.null(parDEC)) {
    #cat("calculating initial values for DEC... \n")
    thetat<- seq(0.1,luDEC,by=.1)
    phit <- seq(0.1,.9,by=.05)
    vect <-merge(phit,thetat,all=T)
    logveroDECv<-function(phitheta){logveroDECs(y, x, z, Nm, Km,time,ind, beta1=beta1, sigmae=sigmae,
    phiDEC=phitheta[1],thetaDEC=phitheta[2], D1=D1,distr=distr, nu=nu,gammas, alphas)}
    logverovec <- apply(vect,1,logveroDECv)
    parDEC <- as.numeric(vect[which.max(logverovec),])
  }
  #phiDEC=parDEC[1]
  #thetaDEC=parDEC[2]
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],parDEC,nu)
  ##
  llji <- logveroDECs(y, x, z,Nm, Km, time,ind, beta1, sigmae,parDEC[1],parDEC[2], D1, distr, nu,gammas,alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")
  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,2*parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","logveroDECs","matrix.sqrt",
                          "dmvnorm","ljtDECs","ljsDECs","ljcnDECs"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovDEC","logveroDECs","matrix.sqrt",
                          "dmvnorm","ljtDECs","ljsDECs","ljcnDECs"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.DEC,objfn = objfn.DEC,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,lb=lb,lu=lu,luDEC=luDEC,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.DEC,objfn = objfn.DEC,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,lb=lb,lu=lu,luDEC=luDEC,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  phiDEC<-EMout$par[(p+2+q2)]
  thetaDEC<-EMout$par[(p+3+q2)]
  if (distr=="sn") {
    nu <- NULL
  }# else nu<-EMout$par[-(1:(p+3+q2))]
  thetav=EMout$par
  gammav<-thetav[-(1:(p+3+q2))]
    if (distr=="st"){
      nu=thetav[p+3+q2+1]
      gammav<-thetav[-(1:(p+3+q2+1))]
      }
  if (distr=="ss"){
      nu=thetav[p+3+q2+1]
      gammav<-thetav[-(1:(p+3+q2+1))]
      }
  if (distr=="scn"){
      nu=thetav[(p+3+q2+1):(p+3+q2+2)]
      gammav<-thetav[-(1:(p+3+q2+2))]
      }


auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }  

  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDECs,y=y, x=x, z=z,Nm=Nm,time=time, beta1=beta1, D1=D1,
                             sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)),  ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y,time=time, x=x, z=z,Nm=Nm, beta1=beta1, D1=D1,
               sigmae=sigmae,phi=c(phiDEC,thetaDEC),depStruct="DEC", distr=distr,nu=nu,gammas=gammas,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiDEC,thetaDEC,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phi1DEC","phi2DEC",
                                   names_dd)
  else names(theta)<- c(colnames(x),"sigma2","phi1DEC","phi2DEC",
                        names_dd,paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=c(phiDEC,thetaDEC),dsqrt=dd,D=D1,gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  MI=InfmatrixDEC(y,x,z,Nm,Km,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,lambda=rep(0,q1),distr = distr,
                           nu = nu,diagD=diagD,skewind=rep(0,q1),gammas,alphas)
   obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)
    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}


DAAREM.CAR1 <- function(formFixed,formFixedNL,formRandom,data,groupVar,timeVar,
                        distr,beta1,sigmae,phiCAR1,D1,nu, gammas,lb,lu,diagD,
                        precisao,informa,max.iter,showiter,showerroriter,
                        algorithm="EM", #"EM" or "DAAREM"
                        parallelphi,parallelnu,ncores,
                        control.daarem=list(),alphas,nknots){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-length(as.numeric(table(ind)))
  ####################################### Clecio #################################
    XN <- model.matrix(formFixedNL, data = data)
    qnl=ncol(XN)
    ni=as.numeric(table(ind))
    #nknots=max(min(ni),3)
    na1=rep(0,m)
    for (i in 1:m){
        aux=ni[1:i]
        na1[i]=sum(aux)
        }
    na2=c(0,na1)
    Nm=list()
    Km=list()
    for (j in 1:qnl){
         Nj=c()
         for (i in 1:m){
              aux=(na2[i]+1):(na2[i+1])
              tij=XN[aux,j]
              Nij=bsplinec(tij,nknots[j],4)
              Nj=rbind(Nj,Nij)              
              }
         Nm[[j]]=Nj
         Dk=diff(diag(ncol(Nj)),differences=2)
         Kj=t(Dk)%*%Dk
         Km[[j]]=Kj
         }
# indg=c()
# for (j in 1:length(Km)){
#      aux=rep(j,ncol(Km[[j]]))
#      indg=c(indg,aux)
#      }
  ng=length(gammas)
  medg=0
  for (i in 1:ng){
       Ni=Nm[[i]]
       mui=Ni%*%gammas[[i]]
       medg=medg+mui
       }         
################################################################################
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(phiCAR1) && length(phiCAR1)!=1) stop("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCAR1)) if (phiCAR1>=1 || phiCAR1<=0) stop ("0<initialValue$phi<1 needed")

  if (is.null(phiCAR1)) {
    lmeCAR = try(lme(formFixed,random=~1|ind,data=data,correlation=corCAR1(form = ~time)),silent=T)
    if (is(lmeCAR,"try-error")) phiDEC =abs(as.numeric(pacf(y-x%*%beta1-medg,lag.max=1,plot=F)$acf))
    else {
      phiDEC = capture.output(lmeCAR$modelStruct$corStruct)[3]
      phiDEC = as.numeric(strsplit(phiDEC, " ")[[1]])
    }
  } else phiDEC <- phiCAR1

  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiDEC,nu)
  ##
  llji <- logveroCAR1s(y, x, z,Nm, Km, time,ind, beta1, sigmae,phiDEC, D1, distr, nu,gammas,alphas)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","logveroCAR1s","matrix.sqrt",
                          "dmvnorm","ljtCAR1s","ljsCAR1s","ljcnCAR1s"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovDEC","logveroCAR1s","matrix.sqrt",
                          "dmvnorm","ljtCAR1s","ljsCAR1s","ljcnCAR1s"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.CAR1,objfn = objfn.CAR1,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.CAR1,objfn = objfn.CAR1,
                    y=y,x=x,z=z,Nm=Nm,Km=Km,time=time,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,diagD=diagD,gammas=gammas,alphas=alphas)
  }
  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  phiDEC<-EMout$par[(p+2+q2)]
  if (distr=="sn") {
    nu <- NULL
  } #else nu<-EMout$par[-(1:(p+2+q2))]
  thetav=EMout$par
  gammav<-thetav[-(1:(p+2+q2))]
    if (distr=="st"){
      nu=thetav[p+2+q2+1]
      gammav<-thetav[-(1:(p+2+q2+1))]
      }
  if (distr=="ss"){
      nu=thetav[p+2+q2+1]
      gammav<-thetav[-(1:(p+2+q2+1))]
      }
  if (distr=="scn"){
      nu=thetav[(p+2+q2+1):(p+2+q2+2)]
      gammav<-thetav[-(1:(p+2+q2+2))]
      }
 
auxg=0
gammas=list()
for (j in 1:qnl){
     dNj=ncol(Nm[[j]])
     gammas[[j]]=gammav[(auxg+1):(auxg+dNj)]
     auxg=auxg+dNj
     }   
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDECs,y=y, x=x, z=z,Nm=Nm,time=time, beta1=beta1, D1=D1,
                             sigmae=sigmae,phiDEC=phiDEC,thetaDEC=1, distr=distr,nu=nu,gammas=gammas,simplify = FALSE)), ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y,time=time, x=x, z=z,Nm=Nm, beta1=beta1, D1=D1,
               sigmae=sigmae,phi=c(phiDEC,1),depStruct="DEC", distr=distr,nu=nu,gammas=gammas,
               simplify = TRUE)
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiDEC,dd,nu)

  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCAR1",names_dd)
  else names(theta)<- c(colnames(x),"sigma2","phiCAR1",names_dd,
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=phiDEC,dsqrt=dd,D=D1,gammas=gammas),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  MI=InfmatrixCAR1s(y,x,z,Nm,Km,time,ind,beta1,sigmae,phiDEC,D1,distr = distr, nu = nu,gammas,alphas)
   obj.out$cov.theta=ginv(MI)
  if (informa) {                 
    desvios<-try(sqrt(diag(ginv(MI))),silent=T)
    if (is(desvios,"try-error")) {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################



