

library(Matrix)

############################### set up scale for a map
mapsetup<- function(xrange,yrange,scale){  #scale is m per inch

  px   <- diff(xrange)/scale  py   <- diff(yrange)/scale  par(pin=c(px,py))

}

############################### map speciesmap_species <- function(specs){

  ns <- length(specs)
  A <- F  for(i in 1:ns){   #extract rows for the ith species   temp <- treedata[treedata[,'species']==specs[i],]   #draw trees for species i   symbols(temp[,'x'],temp[,'y'],circles=temp[,'dbh']/10,inches=F,
    xlim=mapx,ylim=mapy,fg=i,add=A)   #after the 1st species, new species will be added to the existing plot   A <- T  }
}

############################## empirical distribution functionmy.ecdf <- function(x){  n  <- length(x)  xx <- sort(unique(x))                 # find unique values of x  fx <- cumsum(tabulate(match(x,xx))/n) # number of observations for each x  list(x = xx,fx = fx)                  # return the unique values and edcf}
############################# sample plots from a map
sample_plots <- function(mapx,mapy,mapscale,PLOTIT = T){
  yt      <- seq(mapy[1],mapy[2],by=wide)      #y grid locations  yl      <- length(yt)  xt      <- seq(mapx[1],mapx[2],by=wide)      #x grid locations
  xl      <- length(xt)  mapgrid <- cbind(rep(xt,each=yl),rep(yt,xl)) #x and y locations  loc     <- mapgrid[sample(yl*xl,nplot,replace=F),] #samples from grid  sl <- .5*wide                             #plot edges  xbound <- cbind((loc[,1]-sl),(loc[,1]+sl))  ybound <- cbind((loc[,2]-sl),(loc[,2]+sl))  xindex <- c(1,2,2,1,1)  yindex <- c(1,1,2,2,1)  plot.data  <- numeric(0)                   #list of observations  table.spec <- matrix(0,nrow=nspec,ncol=nplot)  rownames(table.spec) <- specname  if(PLOTIT)plot(-1000,0,xlim=mapx,ylim=mapy)  for(i in 1:nplot){  # extract trees on sample plot i and store them in table.spec    xt <- !is.na(cut(treedata[,'x'],breaks=xbound[i,],exclude=NA))    yt <- !is.na(cut(treedata[,'y'],breaks=ybound[i,],exclude=NA))    plot.data      <- append(plot.data,list(treedata[xt & yt,]))    table.spec[,i] <- table(plot.data[[i]][,'species'])   # plot the trees in this plot    if(PLOTIT & dim(plot.data[[i]])[1] > 0){       symbols(plot.data[[i]][,'x'],plot.data[[i]][,'y'],             circles=plot.data[[i]][,'dbh']/20,inches=F,add=T)    }  # draw a box around each plot
    if(PLOTIT){      xvec <- xbound[i,xindex]      yvec <- ybound[i,yindex]      lines(xvec,yvec)
    }  }

  list(specbyplot = table.spec, plotdata = plot.data)
}
######################################

inData <- function(filename, xnames = NULL, ynames = NULL, tname = NULL, 
                   iname = NULL, na.rm = F, INTERCEPT = F){  
                   	
  #read in data file, return design matrix x, response y
  #xnames, ynames, tname, iname are column headings in filename
  #time indicator t
  #individual indicator i

  data <- read.table(filename,header=T)
  
  if(is.atomic(xnames)){
    x    <- data[,xnames]
    if(!is.matrix(x))x <- as.matrix(x)
    if(INTERCEPT){
      intercept <- rep(1,nrow(data))
  	   x <- cbind(intercept,x)
  	   colnames(x) <- c('intercept',xnames)
    }
  }
  if(is.atomic(ynames)){
    y <- data[,ynames]
    if(!is.matrix(y))y <- as.matrix(y) 
    y  <- matrix(y,nrow(data),length(ynames))
    colnames(y) <- ynames
  }
  
  wf <- c(1:nrow(data))
  
  if(na.rm){
  	 wf <- which(is.finite(apply(x,1,sum)) & is.finite(apply(y,1,sum)))
  	 x  <- x[wf,]
  	 y  <- y[wf,]
  }
  
  z  <- list(x = x, y = y)
  
  if(is.atomic(tname))z$t <- data[wf,tname]
  if(is.atomic(iname))z$i <- data[wf,iname]
  z
}
####################################################
simX <- function(n,loX,hiX){                #generate design matrix

  k <- length(loX)
  x <- matrix(1,n,k)
  for(j in 1:k)x[,j] <- runif(n,loX[j],hiX[j])
  x
}
####################################################
simY <- function(x,b,LIKE,r = 1,sigma = 0){     #simulate response

  u <- x%*%b

 # if(length(grep('-',LIKE)) == 1){
 # 	 L <- unlist( strsplit(LIKE,'-') )
 # 	 if(L[1] == 'mvnorm')##########
 # }

  LR   <- paste('r',LIKE,sep='')
  LR <- match.fun(LR)

  if(LIKE == 'pois') u <- exp(u)
  if(LIKE == 'binom')u <- inv.logit(u)

  if(LIKE == 'multinom'){
     zs <- apply(exp(u),1,sum)
     z1 <- 1/(1 + zs)
     zm <- exp(u)/ (1 + zs)
     u  <- cbind(zm,z1)
  }

  if(LIKE == 'pois')    y <- LR(nrow(x),u)
  if(LIKE == 'binom')   y <- LR(nrow(x),1,u)
  if(LIKE == 'norm')    y <- LR(nrow(x),u,sqrt(sigma))
  if(LIKE == 'multinom')y <- myrmultinom(1,u)
  y
}
#########################################
myrmultinom <- function(size,p){  

  #n multinomial r.v. for a n by ncol(p) matrix of probs
  #each row of p is a probability vector

  n     <- nrow(p)
  J     <- ncol(p)
  y     <- matrix(0,n,J)
  sizej <- rep(size,n)
  sumj  <- rep(0,n)
  dpj   <- rep(1,n)
  pj    <- p
  wj    <- c(1:n)

  for(j in 1:(J-1)){
    a     <- round(pj[wj,1],10)
    y[wj,j] <- rbinom(length(wj),sizej[wj],a)
    sumj  <- sumj + y[,j]
    sizej <- size - sumj
    dpj   <- dpj - p[,j]
    pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
    wj    <- which(sumj < size,arr.ind=T) 
  }

  y[,J] <- size - apply(y,1,sum)
  y

}


##############################################################

#like_geom <- function(theta) -sum(dgeom((y-1),theta,log=T)) 

get_lf <- function(bounds,FUN){  #gets the likelihood function

  x <- seq(bounds[1],bounds[2],length=1000)
  q <- x*0
  for(i in 1:1000){q[i] <- FUN(x[i])}
  list(x = x, likelihood = q)
}


treebytime <- function(iindex,tindex,var){ #matrix for var with individuals by time
  vmat <- matrix(NA,max(iindex),max(tindex))
  vmat[cbind(iindex,tindex)] <- var
  vmat
}

like_conedata <- function(b){
  g <- b*x^2
  -sum(dpois(y,g,log=T))
}

like_coneyear <- function(b){
  g     <- b[year]*diamMature^2
  -sum(dpois(coneMature,g,log=T))
}

like_geom <- function(n){
  q <- 1/mean(n)               #ML estimate
  x <- sum(dgeom(n,q,log=T))   #lnL
  list(theta = q, loglik = x)
}


like_weib <- function(param){  #Weibull likelihood  -sum(dweibull(y,param[1],1/param[2],log=T))}

like_exp_binom <- function(rho){  theta <- exp(-rho*20)  -dbinom(6,30,theta,log=T)}

like_bern_logit <- function(pars){
  
  ltheta <- x %*% pars
  theta  <- inv.logit(ltheta)
  s <- y*log(theta) + (1 - y)*log(1 - theta)
  -sum(s)
}


like_norm <- function(par){
  -sum(dnorm(y,par[1],sqrt(par[2]),log=T))
}

surv_prodlim <- function(surv,life){ #product limit estimator of survival
  ctable <- table(surv,life)         #survival by sample date
  yrsum  <- apply(ctable,2,sum)     
  nt     <- n - cumsum(yrsum)        #no. at risk by date
  nt     <- c(n,nt[-nyr])            #1-yr shift
  dt     <- ctable[1,]              #no. died by date
  lamdat <- dt/nt

  s <- rep(1,nyr)                    #survival function
  for(j in 2:nyr)s[j] <- s[j-1]*(1 - lamdat[j])
  list(lamda = lamdat, survprod = s)
}

surv_weib <- function(b){

  lamda <- b[1]
  c     <- b[2]
  beta  <- b[3:length(b)]
  xb    <- x %*% beta
  h0    <- log(c) + c*log(lamda) + (c-1)*log(life[notc]) #only uncensored
  l0    <- -(lamda*life)^c                                 #all individuals
  h     <- h0 + xb[notc]
  s     <- l0*exp(xb)
  s[notc] <- s[notc] + h
  -sum(s)                      #returns -lnL
}

inv.logit <- function(f)(1/(1 + exp(-f)))

logit <- function(f) log(f/(1 - f))


qprob <- function(p){  #evaluate likelihood for pathogen model: multinomial    theta <- p[1]    phi   <- p[2]    s0    <- p[3]    s1    <- p[4]    q00 <- (1 - theta)*(1 - s0) + theta*(1 - phi)*(1 - s1)    q01 <- (1 - theta)*s0 + theta*(1 - phi)*s1    q10 <- theta*phi*(1 - s1)    q11 <- theta*phi*s1    c(q00,q01,q10,q11)}

like_multinom <- function(pars){   q <- qprob(pars)  -sum(nvec*log(q))}

like_pois_probit <- function(par){  a0 <- par[1]  a1 <- par[2]  b  <- par[3] #diameters for zeros  xd     <- xtmp[ytmp == 0]  qd     <- b*xd^2  theta1 <- pnorm(xd,a0,a1)               #probit  pz     <- sum(log((1-theta1) + theta1*dpois(0,qd))) #trees having cones  xd     <- xtmp[ytmp > 0]  yd     <- ytmp[ytmp > 0]  qd     <- b*xd^2  theta2 <- pnorm(xd,a0,a1)  pk     <- sum(log(theta2) + dpois(yd,qd,log=T))  -pz - pk}

#elephant movement example

gk <- function(locnow,q,ss){   dx   <- distance[locnow,]         pvec <- q*exp(-(dx^2)/ss)   pvec/sum(pvec)}

move <- function(locnow,q,ss){   gt  <- gk(locnow,q,ss)      #prob moving to other patches   wk  <- rmultinom(1,1,gt)    #movement vector   c(1:npatch) [wk == 1]       #extract new patch location}update_z <- function(){  #update unknown locations

    for(t in tindex){   #only simulate hidden states      pk   <- gk(z[t-1],qual[(t-1),],sg)*(1-pg[z[t]])      if(t < nt)pk <- pk *gk(z[t+1],qual[(t+1),],sg)      pk   <- pk/sum(pk)      pv   <- rmultinom(1,1,pk)      z[t] <- pseq[pv == 1]    }    z[tindex]
}

update_p <- function(){  #detection probabilities    p1 <- sp1 + uk            p2 <- sp2 + tabulate(z,nbins=npatch) - uk    rbeta(npatch,p1,p2)
}

update_sigma <- function(){  #movement parameter

    props <- tnorm(1,0,Inf,sg,2)   #proposal for sigma
    psum1 <- 0    psum2 <- 0    for(t in 2:(nt-1)){     psum1 <- psum1 + log(gk(z[t],qual[t,],sg)[z[t+1]])     psum2 <- psum2 + log(gk(z[t],qual[t,],props)[z[t+1]])    }    a <- exp(psum2 - psum1)    z <- runif(1,0,1)    if(z < a) sg <- props
    sg
}

#simulation

fitturn <- function(par){ #fit aphid turn data  b0 <- par[1]  b1 <- par[2]  b2 <- par[3]  mu <- (b0 + b1*cos(b2 + x))/2/pi/b0  #mean function  -sum(dpois(y,2*pi/m*mu,log=T))       # - ln likelihood}

rturn <- function(n,c,cf){  #generate random deviates  x  <- numeric(0)  nx <- length(x)  while(nx < n){    q  <- runif(n-nx,-pi,pi)    a  <- (b0 + b1*cos(b2 + q))/c/cf                #unnormalized    z  <- runif(n,0,1)    xz <- q[z < a]    x  <- c(xz,q)    nx <- length(x)  }  x}

rtrunc <- function(n,lo,hi,p1,p2,FN){    #truncated using inv dist sampling   #truncated 2 parameter distribution	  pf <- paste('p',FN,sep='')  qf <- paste('q',FN,sep='')  pf1 <- match.fun(pf)  qf1 <- match.fun(qf)  z  <- runif(n,pf1(lo,p1,p2),pf1(hi,p1,p2))  qf1(z,p1,p2)}


tnorm <- function(n,lo,hi,mu,sig){   

  #normal truncated lo and hi

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig)

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}


acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H	# if BLOCK, then accept as a block,	# otherwise, accept individually  nz          <- length(x0)  #no. to accept  if(BLOCK)nz <- 1    a    <- exp(p1 - p0)       #acceptance PR  z    <- runif(nz,0,1)  keep <- which(z < a,arr.ind=T)    if(BLOCK & length(keep) > 0)x0 <- x1  if(!BLOCK)                  x0[keep] <- x1[keep]             ac <- length(keep)          list(x = x0, accept = ac)}

pmake <- function(pars){    p   <- pars[1]                   #frequency of a  f   <- pars[2]                   #inbreeding coefficient   paa <- p^2 + f*p*(1 - p)	  #frequency of p.aa  pab <- 2*p*(1 - p)*(1 - f)	  #frequency of p.ab  pbb <- (1 - p)^2 + f*p*(1 - p)  #frequency of p.bb  c(paa,pab,pbb)}

minf <- function(p){   #minimum inbreeding coefficient  lof <- rep(-1,length(p))  lop <- -(1 - p)/p  hip <- -p/(1 - p)  lof[p > (1 - p)] <- lop[p > (1 - p)]  lof[p < (1 - p)] <- hip[p < (1 - p)]  lof}pf.update <- function(p,f){  #log posterior  multlik(c(p,f)) + dbeta(p,g1,g2,log=T) +            dnorm(f,priorF,priorFSD,log=T)}
multlik <- function(pars){  #multinom likelihood   p    <- pmake(pars)  pmat <- matrix(rep(p,npop),3,npop)  sum(y*log(pmat))}


tnorm.mvt <- function(avec,muvec,smat,lo,hi){   

  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 

 # smat <- nearPD(smat)$mat

 # avec    <- muvec
  avec[1] <- tnorm(1,lo[1],hi[1],muvec[1],sqrt(smat[1,1]))

  for(k in 2:length(muvec)){
    skk <- smat[-k,-k]
  #  skk <- as.matrix(nearPD(smat[-k,-k])$mat)

     testv <- try(chol(skk),T)
     if(inherits(testv,'try-error')){
      
         avec[k] <- 0
         next

      }

    piece1 <- smat[-k,k] %*% chol2inv(testv)
    muk <- muvec[k] + piece1 %*% (avec[-k] - muvec[-k])
    sgk <- as.numeric(smat[k,k] - piece1 %*% smat[k,-k])
    if(sgk < .000000001)sgk <- .000000001
    avec[k] <- tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
  }
  avec
}


update_pf <- function(){

  propp  <- tnorm(1,.02,.98,pg,.002) #propose pa  fl     <- minf(propp)              #lower f limit for pa  propf  <- tnorm(1,fl,1,fg,.01)     #propose f > fl  pnew   <- pf.update(propp,propf)  pnow   <- pf.update(pg,fg)  atmp   <- acceptMH(pnow,pnew,c(pg,fg),c(propp,propf),BLOCK=T)  pg     <- atmp$x[1]  fg     <- atmp$x[2]  ac     <- atmp$accept

  list(pg = pg, fg = fg, accept = ac)
}

##sample from posteriors for a linear regression

b.update <- function(){  V <- solve(crossprod(x)/sg + priorIV)  v <- crossprod(x,y)/sg + priorIV %*% priorB  t(rmvnorm(1,V%*%v,V))}

v.update <- function(){  u1 <- s1 + n/2  u2 <- s2 + .5*crossprod(y - x%*%bg)  1/rgamma(1,u1,u2)}

#######################################################

gibbsLoop <- function(LIKE,ng,x,y,b,sigma = 0){

  k      <- length(b)
  bgibbs <- matrix(NA,ng,k);  colnames(bgibbs) <- paste('b',c(1:k),sep='-')
  sgibbs <- rep(NA,ng)
  if(sigma == 0)sgibbs <- numeric(0)
  r <- 1
  if(is.matrix(y))r <- ncol(y)
  pred   <- rep(0,nrow(x)*r)
  pred2  <- pred

  bg <- b
  sg <- sigma

  for(g in 1:ng){

    bg <- bUpdateGibbs(x,y,bg,LIKE,sg)
    if(sigma > 0){
    	 sg <- sigmaUpdate(x,y,bg)
    	 sgibbs[g]  <- sg
    }

    py    <- as.vector(simY(x,bg,LIKE,r,sg))
    pred  <- pred + py
    pred2 <- pred2 + py^2

    bgibbs[g,] <- bg
  }

  ymean <- pred/ng
  yse   <- sqrt(pred2/ng - ymean^2)

  list(bgibbs = bgibbs,sgibbs = sgibbs, ymean = ymean, yse = yse)
}
####################################################
bUpdateGibbs <- function(x,y,b,LIKE,sigma = 0){

  if(LIKE == 'norm')    return( bUpdateNorm(x,y,sigma) )
  if(LIKE == 'multinom')return( bUpdateMNom(x,y,b) )

  b <- matrix(b,length(b),1)
  c <- t(myrmvnorm(1,t(b),pBVar))

  znow <- x%*%b
  znew <- x%*%c

  if(LIKE == 'pois'){
     pnow <- dpois(y,exp(znow),log=T)
     pnew <- dpois(y,exp(znew),log=T)
  }
  if(LIKE == 'binom'){
     pnow <- dbinom(y,1,inv.logit(znow),log=T)
     pnew <- dbinom(y,1,inv.logit(znew),log=T)
  }

  pnow <- sum(pnow) + mydmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(t(c),priorB,priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}
####################################################
bUpdateMNom <- function(x,y,b){

  bvec <- as.vector(b)

  tmp  <- bmultiProp(b)
  c    <- tmp$c
  cvec <- tmp$cvec

  z <- x%*%b
     zs    <- apply(exp(z),1,sum)
     z1    <- 1/(1 + zs)
     zm    <- exp(z)/ (1 + zs)
     znow  <- cbind(zm,z1)

  z <- x%*%c
     zs   <- apply(exp(z),1,sum)
     z1   <- 1/(1 + zs)
     zm   <- exp(z)/ (1 + zs)
     znew <- cbind(zm,z1)

  pnow <- sum(y*log(znow)) + mydmvnorm(bvec,priorB,priorVB,log=T)
  pnew <- sum(y*log(znew)) + mydmvnorm(cvec,priorB,priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}

####################################################
bmultiProp <- function(b = matrix(0,k,r-1)){  

    bvec <- as.vector(b)
    cvec <- myrmvnorm(1,t(bvec),pBVar)
    c    <- matrix(cvec,nrow(b),ncol(b))

  list(c = c, cvec = cvec)
}

####################################################
bUpdateNorm <- function(x,y,sigma){

  V <- solve(crossprod(x)/sigma + priorIV)
  v <- crossprod(x,y)/sigma + priorIV %*% priorB
  t(myrmvnorm(1,t(V%*%v),V))
}
####################################################
sigmaUpdate <- function(x,y,b){

  u1 <- s1 + nrow(x)/2
  u2 <- s2 + .5*sum( (y - x%*%b)^2 )
  1/rgamma(1,u1,u2)
}
##########################################################
processPars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                        sigOnly = F,burnin=1,xlimits = NULL){  

  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb[,-wi]
    if(length(xtrue) > 0)xtrue <- xtrue[-wi]

    wq   <- apply(btmp,2,quantile,c(.025,.975))  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    xgb  <- btmp[,-wq]
    if(length(xtrue) > 0)xtrue <- xtrue[-wq]
   }

  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1)     xgb <- xgb[-c(1:burnin),]
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)

  out <- t(rbind(apply(xgb,2,mean),apply(xgb,2,quantile,c(.025,.975))))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('mean','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('mean','0.025','0.975','true value')
  }

  armat <- matrix(0,nc,10)  #for AR model
  
  
 # for(j in 1:nc)armat[j,] <- ar(xgb[,j],aic=F,order.max = 10)$ar[1:10]

 # if(!is.null(colnames(xgb)))rownames(armat) <- colnames(xgb)
#  colnames(armat) <- paste('AR',c(1:10),sep='-')

  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc))

  if(CPLOT){
      for(j in 1:nc){
       plot(xgb[,j],type='l')
       abline(h=out[j,],lty=2)
       if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
       title(colnames(xgb)[j])
     }
  }
  xlims <- xlimits
  if(DPLOT){
      for(j in 1:nc){
        xj <- density(xgb[,j])
        if(is.null(xlimits))xlims <- range(xj$x)
        plot(xj$x,xj$y,type='l',xlim=xlims)
        abline(v=out[j,],lty=2)
        if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
        title(colnames(xgb)[j])
     }
  }
  list(summary = signif(out,4), ar_model = signif(armat,3))

}


gibbsStates <- function(time,x,y,xtrue=numeric(0)){  nc <- ncol(x)  par(mfrow=c(1,1))  plot(time,y)  xci <- apply(x,2,quantile,c(.5,.025,.975))  lines(time,xci[1,])  lines(time,xci[2,],lty=2)  lines(time,xci[3,],lty=2)  if(length(xtrue) > 0)lines(time,xtrue,col='red')}
#############state space models

SSRWupdateX <- function(){         #update x's, random walk


  for(t in 1:nt){

    VI <- 0
    v  <- 0

    if(!t %in% wm){          #observations
      v  <- y[t]/tg
      VI <- 1/tg
    }

    if(t < nt){              #t+1 term excluded for last 
      v  <- v + xg[t+1]/sg 
      VI <- VI + 1/sg
    }

   if(t > 1){                #t-1 term excluded for 1st 
      v  <- v + xg[t-1]/sg
      VI <- VI + 1/sg
   }

   V     <- 1/VI
   xg[t] <- rnorm(1,V*v,sqrt(V))
  }
  xg
}

SSupdateS <- function(FUN,...){        #process error
  
  FUN <- match.fun(FUN)
  x2 <- xg[-1]
  x1 <- xg[-nt] +  FUN(...)
  u1 <- s1 + (nt-1)/2
  u2 <- s2 + .5*sum( ((x2 - x1)^2) )
  1/rgamma(1,u1,u2)
}

SSupdateT <- function(){

  u1 <- v1 + (nt-nm)/2
  u2 <- v2 + .5*sum( (y[-wm] - xg[-wm])^2)
  1/rgamma(1,u1,u2)
}

SSRWupdateB <- function(){  #regression parameters in rw model  V <- 1/( (nt - 1)/sg + 1/vb )  v <- (x[nt] - x[1])/sg + b0/vb  rnorm(1,V*v,sqrt(V))}
SSLGupdateB <- function(){  z <- xg[-1] - xg[-nt]  X <- cbind(rep(1,(nt-1)),exp(xg[-nt]))  V <- solve(crossprod(X)	+priorIVb)  v <- crossprod(X,z) + priorIVb%*%priorb  rmvnorm(1,V%*%v,V)  }
SSMetUpdateX <- function(){     #Metropolis for state-space  for(t in 1:nt){    xprop <- rnorm(1,xg[t],.1)    pnow <- 0    pnew <- 0    if(!t %in% wm){          #observations      pnow  <- pnow + dnorm(y[t],xg[t],sqrt(tg),log=T)      pnew  <- pnew + dnorm(y[t],xprop,sqrt(tg),log=T)    }    if(t < nt){              #t+1 term excluded for last       pnow  <- pnow + dnorm(xg[t+1],fx(xg[t]),sqrt(sg),log=T)      pnew  <- pnew + dnorm(xg[t+1],fx(xprop),sqrt(sg),log=T)    }   if(t > 1){                #t-1 term excluded for 1st       pnow  <- pnow + dnorm(xg[t],fx(xg[t-1]),sqrt(sg),log=T)      pnew  <- pnew + dnorm(xprop,fx(xg[t-1]),sqrt(sg),log=T)   }   xg[t] <- acceptMH(pnow,pnew,xg[t],xprop,BLOCK=T)$x   }  xg}

SSLGupdateRK <- function(){  prop <- tnorm(2,c(0,0),c(10,1000),c(rg,kg),c(.01,1))   rp   <- prop[1]; kp <- prop[2]  nnow <- xg[-nt] + fx(xg[-nt],c(rg,kg,dt))   #predicted n  nnew <- xg[-nt] + fx(xg[-nt],c(prop,dt))	  pnow <- sum(dnorm(xg[-1],nnow,sqrt(sg*dt),log=T)) +          sum(dnorm(c(rg,kg),c(mur,muk),c(sdr,sdk),log=T))  pnew <- sum(dnorm(xg[-1],nnew,sqrt(sg*dt),log=T)) +	       sum(dnorm(c(rp,kp),c(mur,muk),c(sdr,sdk),log=T))  atmp <- acceptMH(pnow,pnew,c(rg,kg),prop,BLOCK=T)  rg   <- atmp$x[1]  kg   <- atmp$x[2]  list(rg = rg, kg = kg, aa = atmp$accept)}SSLGupdateX <- function(){   #propose/accept x as a block  xprop <- tnorm(nt,0,1000,xg,.2)                     #proposed x  xpnew <- xprop[-nt] + fx(xprop[-nt],c(rg,kg,dt))    #mean proposed  xpnow <- xg[-nt]    + fx(xg[-nt],c(rg,kg,dt))       #mean current  pnow <- sum(dnorm(xg[-1],xpnow,sqrt(sg*dt),log=T)) +	       sum(dnorm(y[-wm],xg[-wm],sqrt(tg),log=T))  pnew <- sum(dnorm(xprop[-1],xpnew,sqrt(sg*dt),log=T)) +	       sum(dnorm(y[-wm],xprop[-wm],sqrt(tg),log=T))  atmp <- acceptMH(pnow,pnew,xg,xprop,BLOCK=T)  xg   <- atmp$x  list(xg = xg, aa = atmp$accept)}SSLGupdateS <- function(){    xp <- xg[-nt] + fx(xg[-nt],c(rg,kg))*dt    #mean current  u1 <- s1 + (nt-1)/2  u2 <- s2 + .5*sum( ((xg[-1] - xp)^2)/dt )  1/rgamma(1,u1,u2)}myrmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean))){
    ev <- eigen(sigma, sym = TRUE)$values
    if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1])))
        warning("sigma is numerically not positive definite")
    sigsvd <- svd(sigma)
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval + mean
}
mydmvnorm <- function(x,mean,sigma,log=FALSE){

  #mv normal density

    if (is.vector(x))x <- matrix(x, ncol = length(x))
    if (is.vector(mean))mean <- matrix(mean, ncol = length(x))

    distval <- mahalanobis(x, mean, sigma)
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    if(log)return(logretval)
    exp(logretval)
}

ZIPparUpdate <- function(){

  cprop <- t(myrmvnorm(1,c(ag,bg),parcov))
  ap <- cprop[1:nv]
  bp <- cprop[(nv+1):(nv+nx)]

  if(DIST == 'Poisson'){
    pnow  <- -ziPoisNegLik(c(ag,bg))
    pnew  <- -ziPoisNegLik(c(ap,bp))
  }
  if(DIST == 'LN'){
    pnow  <- -ziLNNegLik(c(ag,bg,sg))
    pnew  <- -ziLNNegLik(c(ap,bp,sg))
  }

  pnow  <- pnow + mydmvnorm(t(ag),aprior,aVar,log=T) +
                  mydmvnorm(t(bg),bprior,bVar,log=T)
  pnew  <- pnew + mydmvnorm(t(ap),aprior,aVar,log=T) +
                  mydmvnorm(t(bp),bprior,bVar,log=T)
  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a){
    ag <- ap
    bg <- bp
  }

  ss <- 0
  if(DIST == 'LN'){
    ss <- sg
    sp <- tnorm(1,0,10,sg,.005)
    pnow  <- -ziLNNegLik(c(ag,bg,sg))
    pnew  <- -ziLNNegLik(c(ag,bg,sp))
    a <- exp(pnew - pnow)
    z <- runif(1,0,1)
    if(z < a)ss <- sp
  }
     
  list(ag = ag, bg = bg,sg = ss)
}

ziPoisNegLik <- function(par){

  a <- matrix(par[1:nv],nv,1)
  b <- matrix(par[(nv+1):(nv+nx)],nx,1)

  VA   <- vb%*%a
  LB   <- xb%*%b

  theta <- inv.logit(VA)
  lamda <- exp(LB)

  pnow  <- rep(0,n)
  pnow[yb == 0] <- log(1 - theta + theta*dpois(0,A*lamda))[yb == 0]
  pnow[yb > 0]  <- log(theta[yb > 0]) + dpois(yb[yb > 0],(A*lamda)[yb > 0],log=T)
  -sum(pnow)
}
#########################################
ziLNNegLik <- function(par){

  a <- matrix(par[1:nv],nv,1)
  b <- matrix(par[(nv+1):(nv+nx)],nx,1)
  s <- par[length(par)]

  VA   <- vb%*%a
  LB   <- xb%*%b

  theta <- inv.logit(VA)

  pnow  <- rep(0,n)
  pnow[yb == 0] <- log(1 - theta)[yb == 0]
  pnow[yb > 0]  <- log(theta[yb > 0]) + dnorm(log(yb[yb > 0]),LB[yb > 0],s,log=T)
  -sum(pnow)
}

####################################################
bgsampMVN <- function(){  #sample reg pars for MVN 
	
  V   <- solve(crossprod(x)) 
  v   <- crossprod(x,y)
  mu  <- V%*%v
  vaa <- kronecker(sg,V)   
  
 # vaa <- nearPD(vaa)$mat
  matrix( myrmvnorm(1,as.vector(mu),vaa) ,p,ns,byrow=F)

}
######################################
wishsamp <- function(){   #sample from Inv Wishart

   scp  <- crossprod((y - x %*% bg))
   vmat <- solve(scp + prior.W*prior.WDF)
   v2   <- ns + prior.WDF
   stmp <- myrmvnorm(v2,matrix(0,v2,ns),vmat)
   crossprod(stmp)

}
####################################################
ysampMVNPois <- function(){  #sample y's for MVN 1st stage, Pois 2nd stage

  propy <- matrix(rnorm(length(y),y,.01),n,ns,byrow=F)
  
  pnow <- rep(0,n)
  pnew <- rep(0,n)
  
  for(i in 1:n){
    pnow[i] <- mydmvnorm(y[i,],(x %*% bg)[i,],sg,log=T)
    pnew[i] <- mydmvnorm(propy[i,],(x %*% bg)[i,],sg,log=T)
  }
  
  pnow <- pnow + rowSums(dpois(z,Ag*exp(y),log=T))
  pnew <- pnew + rowSums(dpois(z,Ag*exp(propy),log=T))

  a <- exp(sum(pnew) - sum(pnow))
  zz <- runif(1,0,1)
  accept <- 0
  if(zz < a){
     y <- propy
     accept <- 1
  }
  list(y = y, a = accept)
  
}
#########################

sampleEffort <- function(){  #if Poisson effort is unknown
	
	u1 <- a1 + apply(z,1,sum)
	u2 <- a2 + apply(exp(y),1,sum)
	a <- rgamma(n,u1,u2)
	matrix(a,n,ns)
}

#########################
simMVNData <- function(){
	
   covars <- c('intercept',paste('x',c(1:(p-1)),sep='-') )
   p      <- length(covars)
   x      <- matrix(runif(n*p,-1,1),n,p)
   x[,1]  <- 1
   colnames(x) <- covars
   
   b           <- matrix(runif(p*ns,-2,1),p,ns)
   specnames   <- paste('S',c(1:ns),sep='-')
   colnames(b) <- specnames
   rownames(b) <- covars
   sigma <- diag(runif(ns,0,.1),ns)
   rownames(sigma) <- specnames
   colnames(sigma) <- specnames

   mu <- x %*% b
   y     <- myrmvnorm(n,mu,sigma)
   sigma <- cov(y)
   y     <- myrmvnorm(n,mu,sigma)
   colnames(y) <- specnames
   
   list(x = x, b = b, y = y, s = sigma, specnames = specnames)
}
YtoZ <- function(y){ 
	#multivar logit to fractions
  zs   <- apply(exp(y),1,sum)
  z1   <- 1/(1 + zs)
  zm   <- exp(y)/ (1 + zs)
  cbind(zm,z1)
}

ZtoY <- function(z){     
	#fractions to multivar logit
 # log(z[,-(ns+1)]*(1 + z))       
  
  log(z[,-(ns+1)]/
  (1 - apply(z[,-(ns+1)],1,sum)))
  
}

predYMVN <- function(x,bg,sg){ 
	#predict multivariate y

  vaa <- kronecker(diag(1,n),sg)
  matrix(myrmvnorm(1,as.vector(x%*%bg),vaa),n,ns)
}

#####################
distmat <- function(xt,yt,xs,ys){
    xd <- outer(xt,xs,function(xt,xs) (xt - xs)^2)
    yd <- outer(yt,ys,function(yt,ys) (yt - ys)^2)
    t(sqrt(xd + yd)) 
}
#####################
getSetup <- function(mapx,mapy,grid){   #setup map for movement

  nx   <- diff(mapx)/grid + 1
  ny   <- diff(mapy)/grid + 1
  xs <- seq(mapx[1],mapx[2],length=nx)
  ys <- seq(mapy[1],mapy[2],length=ny)

  pj   <- as.matrix(expand.grid(x=seq(0,100,length=nx),y=seq(0,100,length=ny)))
  dj   <- distmat(pj[,1],pj[,2],pj[,1],pj[,2])

  list(xs = xs, ys =  ys, pj = pj, dj = dj)
}
#################################
getNeighbors <- function(j){   #find neighbor cells for movement model

  dz     <- matrix(distj[j,],nrow=length(j),byrow=F)
  nb     <- t(apply(dz,1,order,decreasing=F))[,1:nm]
  dindex <- cbind(rep(j,each=nm),as.vector(t(nb)))
  distNb <- matrix(distj[dindex],nrow=length(j),byrow=T)

  list(nb = nb, di = distNb)
}
##########################

MoveUpdateB <- function(){  #update reg pars for movement model

   bp   <- t(myrmvnorm(1,t(bg),parcov))

   pnow <- 0
   pnew <- 0

   for(t in 1:(nt-1)){

     znext <- which(hood[t,] == z[t+1])
     xt    <- cbind(xvar[hood[t,]] - xvar[z[t]],distn[t,]) #subtract mean to get gradient
     e     <- exp(xt%*%bg)
     th    <- (e/sum(e))[znext]
     pnow  <- pnow + log(th)

     e     <- exp(xt%*%bp)
     th    <- (e/sum(e))[znext]
     pnew  <- pnew + log(th)

  }

  pnow <- pnow + mydmvnorm(t(bg),bprior,bVar,log=T)
  pnew <- pnew + mydmvnorm(t(bp),bprior,bVar,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)bg <- bp
  bg
}



