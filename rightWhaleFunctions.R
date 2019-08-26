healthIndex <- function(CNAME){  #indices from 3 (healthy) to 1 (severe problems)

  if(!CNAME %in% c('BodyFatCode', 'LeftRakeMarkCode', 'RightRakeMarkCode',
                   'SkinCode', 'CyamidCode')) stop ('Covariate not recognized')

  dataBF  <- sights[,CNAME]
  dataBF[dataBF == 'X'] <- NA
  dataBF  <- as.numeric(as.character(dataBF))

  dtmp <- dataBF

  if(CNAME %in% c('BodyFatCode','LeftRakeMarkCode','RightRakeMarkCode')){
    dtmp[dataBF == 3] <- 1
    dtmp[dataBF == 1] <- 3
    dataBF <- dtmp
  }else{
    dtmp[dataBF == 2] <- 1
    dtmp[dataBF == 1] <- 2
    dataBF <- dtmp
  }

  dataBF
}



updateLamda <- function(){

  a1 <- rowSums(sightMat)
  b1 <- rowSums(effortState*rState)

  a <- aLam + a1
  b <- bLam + b1
  rgamma(n,a,b)
}

updateLamda <- cmpfun(updateLamda)


updateAlam <- function(){

  anew <- tnorm(1,0,10,aLam,.006)

  pnb <- effortState*rState/(effortState*rState + bLam)

  wp <- which(pnb > 0,arr.ind=T)

  py <- dnbinom(sightMat[wp],size=aLam,pnb[wp],log=T)
  py[is.na(py)] <- 0
  pnow <- sum(py) + dgamma(aLam,aLam1,bLam1,log=T)

  py <- dnbinom(sightMat[wp],size=anew,pnb[wp],log=T)
  py[is.na(py)] <- 0
  pnew <- sum(py) + dgamma(anew,aLam1,bLam1,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)aLam <- anew
  aLam
}
updateAlam <- cmpfun(updateAlam)


getMProb <- function(sex,w1,w2,month){
  pm <- femMoveProb[w1,w2,monNow]
  if(sex == 'M')pm <- maleMoveProb[w1,w2,monNow]
  if(sex == 'X')pm <- xMoveProb[w1,w2,monNow]
  if(length(w1) > 1 & length(w2) > 1){
    wd <- which(dim(pm) == nreg)
    pm <- apply(pm,wd,sum) - apply(pm,wd,prod)
  }
  pm
}


getMProbMat <- function(sex,w1,w2,month, juv = rep(0, length(w1))){

  tiny <- 1e-10

  fm <- femMoveProb[,,month]
  mm <- maleMoveProb[,,month]
  jm <- juvMoveProb[,,month]

  fm[fm == 0] <- tiny
  mm[mm == 0] <- tiny
  jm[jm == 0] <- tiny

  if(is.matrix(w2)){              # from vector w1 to matrix w2 of (0,1)
    pm    <- w2*0
    wfem  <- which(sex == 'F' & juv == 0 & w1 > 0)
    wmale <- which(sex == 'M' & juv == 0 & w1 > 0)
    wjuv <- which(juv == 1 & w1 > 0)
    pm[wfem,]  <- t(fm[,w1[wfem]])*w2[wfem,]
    pm[wmale,] <- t(mm[,w1[wmale]])*w2[wmale,]
    pm[wjuv,]  <- t(jm[,w1[wjuv]])*w2[wjuv,]

    return(pm)
  }
  if(is.matrix(w1)){            # from matrix w1 to vector w2
    pm    <- w1*0
    wfem  <- which(sex == 'F' & juv == 0 & w2 > 0)
    wmale <- which(sex == 'M' & juv == 0 & w2 > 0)
    wjuv <- which(juv == 1 & w2 > 0)

    pm[wfem,]  <- fm[w2[wfem],]*w1[wfem,]
    pm[wmale,] <- mm[w2[wmale],]*w1[wmale,]
    pm[wjuv,]  <- jm[w2[wjuv],]*w1[wjuv,]

    return(pm)
  }
  pm <- w1*0                   #from vector to vector
  wfem  <- which(sex == 'F' & juv == 0 & w2 > 0)
  wmale <- which(sex == 'M' & juv == 0 & w2 > 0)
  wjuv <- which(juv == 1 & w2 > 0)
  pm[wfem]  <- fm[cbind(w2[wfem],w1[wfem])]
  pm[wmale] <- mm[cbind(w2[wmale],w1[wmale])]
  pm[wjuv]  <- jm[cbind(w2[wjuv],w1[wjuv])]
  pm

}
getMProbMat <- cmpfun(getMProbMat)


moveSum <- function(from,to,sex,juv){

  wf   <- which(sex == 'F' & juv == 0 & from > 0 & to > 0)
  ftab <- table(regID[to[wf]],regID[from[wf]])

  wm   <- which(sex == 'M' & juv == 0 & from > 0 & to > 0)
  mtab <- table(regID[to[wm]],regID[from[wm]])

  wx   <- which(sex == 'X' & juv == 0 & from > 0 & to > 0)
  xtab <- table(regID[to[wx]],regID[from[wx]])

  wx   <- which(juv == 1 & from > 0 & to > 0)
  jtab <- table(regID[to[wx]],regID[from[wx]])

  list(fem = ftab, male = mtab, x = xtab, juv = jtab)
}
moveSum <- cmpfun(moveSum)


fillSites <- function(www, rt, first, last, t){

  rsum <- rowSums(rt)
  w0   <- which(rsum == 0 & first > 0 & last > 0)
  if(length(w0) == 0) return(rt)

  q1 <- which(myrmultinom(1, matrix(1 / nreg, nrow = length(w0), ncol = nreg)) == 1, arr.ind = T)
  qmat <- matrix(0, length(w0), nreg)
  qmat[q1] <- 1
  rt[w0, ] <- qmat

  ws <- which(sightMat[www, t, ] > 0, arr.ind = T)
  if(length(ws) > 0){
    if(length(www) == 1) ws <- cbind(1, ws)
    s1 <- unique(ws[, 1])
    rt[s1, ] <- 0
    rt[ws]  <- 1
  }

  rt
}
fillSites <- cmpfun(fillSites)


moveSurv <- function(ws0,ws1,ws2,ws3,surv0,surv1,state0,state1,state2,m0,m1,notSight){

  ns    <- ncol(m0)
  plive <- pdead0 <- pdead1 <- m0 * 0

  plive      <- m0 * surv0 * notSight * m1 * surv1          #survived both in each region

  if(length(ws2) > 0){                              #dead at t or at t+1
    plive[ws2,]  <- 0
    pdead1[ws2,] <- m0[ws2,] * surv0[ws2,] * notSight[ws2,] * (1 - surv1[ws2,]) #died on t, t+1
    pdead0[ws2,] <- m0[ws2,] * (1 - surv0[ws2,])                            #died on t-1, t
  }
  if(length(ws3) > 0){                             #known alive after t
    pdead0[ws3,] <- 0
  }
  if(length(ws0) > 0){                             #obs dead at t
    plive[ws0,]  <- 0
    pdead0[ws0,] <- 0
    pdead0[ws0, state1[ws0]] <- 1
    pdead1[ws0,] <- 0
  }
  if(length(ws1) > 0){                             #obs dead at t+1
    pdead0[ws1,] <- 0
    pdead1[ws1,] <- 0
    pdead1[ws1, state1[ws1]] <- 1
  }

  pall <- cbind(plive, pdead1, pdead0)  #(survive both), (survive,died), (died)
  pall <- pall/matrix(rowSums(pall), nrow(pall), ncol(pall))
  pall[which(is.na(pall))] <- 0

  q0 <- myrmultinom(1, matrix(pall,nrow = nrow(plive)))
  q1 <- which(q0 == 1, arr.ind = T)

  q2 <- q1[,1] * 0
  q2[q1[,1]] <- q1[,2]
  dead1 <- which(q2 > ns & q2 <= (2 * ns))
  dead0 <- which(q2 > (2 * ns))

  surv <- rep(1, nrow(m0))
  surv[dead0] <- 0

  stateNew <- q2
  wws <- which(is.finite(stateNew) & stateNew > ns)
  if(length(wws) > 0)stateNew[stateNew > ns] <- stateNew[stateNew > ns] - ns

  wws <- which(is.finite(stateNew) & stateNew > ns)
  if(length(wws) > 0)stateNew[stateNew > ns] <- stateNew[stateNew > ns] - ns

  list(surv = surv, states = stateNew)
}
moveSurv <- cmpfun(moveSurv)


getRS <- function(www, rr, ss, t){   #wts is which individuals sighted at t

  rt <- rr
  if(length(www) == 1)rt <- matrix(rt,1)

  rt <- fillSites(www, rt, t + 1 - firstSight[www], death[www] - t + 1, t)

  wdead <- which(deadID %in% www & deadTime < t)
  if(length(wdead) > 0){
    rt[match(deadID[wdead], www), ] <- 0
  }

  wss <- which(ss > 0)


  rt[wss, ] <- 0
  if(length(wss) == 1) rt[wss, ss[wss]] <- 1
  if(length(wss) > 1) rt[cbind(wss, ss[wss])] <- 1

  rt
}
getRS <- cmpfun(getRS)


getSight <- function(ss){

  if(length(dim(ss)) == 2){
    sout <- matrix(0, 1, nrow(ss))
    wss  <- apply(ss, 1, which.max)
    sss  <- rowSums(ss)
    sout[sss > 0] <- wss[sss > 0]
    return(sout)
  }

  sout <- matrix(0, dim(ss)[1], dim(ss)[2])

  for(tt in c(1:dim(ss)[2])){
    wss <- unique(which(ss[, tt, ] > 0, arr.ind = T)[, 1])
    if(length(wss) == 0)next
    if(length(wss) == 1)wmm <- which.max(ss[wss, tt, ])
    if(length(wss) > 1) wmm <- apply(ss[wss, tt, ], 1, which.max)
    #      if(length(wss) > 1) wmm <- max.col(ss[wss, tt, ])
    sout[wss, tt] <- wmm
  }
  sout
}
getSight <- cmpfun(getSight)

updateStates <- function(){

  #deadID - individual that died
  #deadTime <- when died
  #deadReg <- regions where death sighted

  t1 <- min(firstSight) + 1

  rs <- rState
  rs[sightMat > 0] <- 1

  death[death < lastSight] <- lastSight[death < lastSight]
  if(Survhack){
    death[death > (lastSight + missingInt)] <- lastSight[death > (lastSight + missingInt)] + missingInt
    death[death > (nt + 1)] <- nt + bnt
  }

  if(!subAn){death[deadID] <- deadTime}

  moveMale <- moveFem <- moveX <- moveJ <- maleMove*0

  if(!REGINTERCEPT){
    bnew <- tnorm.mvt(bg[1:2], bg[1:2], propBS[1:2, 1:2], loBS[1:2], hiBS[1:2])
    bnew <- c(bnew[1], rep(bnew[2], nreg))
  }

  if(REGINTERCEPT) {
    bnew <- tnorm.mvt(bg, bg, propBS, loBS, hiBS)
    ssss <- sample(c(1:(nreg + 1)), sample(c(1:(nreg + 1)), 1) )
    bnew[ssss] <- bg[ssss]
  }



  psnow <- psnew <- 0
  #for(i in 1:100){
  for(t in t1:nt){

    www    <- which(firstSight <= t & death >= t)
    nww    <- length(www)
    if(nww == 0)next

    age <- (t - firstSight[www])/12
    ja  <- age*0
    ja[age <= juvAge] <- 1

    monNow <- monYr[t,1]

    if(t < nt) sightT <- getSight(sightMat[www, (t - 1):(t + 1), ])
    if(t == nt)sightT <- cbind(getSight(sightMat[www, (t - 1):t, ]), rep(0, nww))

    effort <- effortmat[t,]

    rt0 <- getRS(www, rs[www, t - 1, ], sightT[, 1], t) # why is this one always get set to state at t if not sighted at t-1, but sighted at t?
    rt1 <- getRS(www, rs[www, t, ], sightT[, 2], t)

    if(t < nt) rt2 <- getRS(www, rs[www, t + 1, ], sightT[, 3], t)
    if(t == nt)rt2 <- getRS(www, rs[www, t, ], sightT[, 2], t)

    rs[www,t,] <- 0

    ws0 <- match(deadID[deadTime == t], www)
    ws1 <- match(deadID[deadTime == (t + 1)], www)
    ws2 <- which(death[www] == t | death[www] == (t + 1))
    ws3 <- which(lastSight[www] > t)

    state0     <- rowMaxs(rt0 * iByr[www, ])                        #location t-1, t, t+1
    state1     <- rowMaxs(rt1 * iByr[www, ])
    state2     <- rowMaxs(rt2 * iByr[www, ])

    rnow       <- rt1 * 0 + 1
    ww2        <- which(sightT[, 2] > 0)
    rnow[ww2, ]  <- 0
    rnow[cbind(ww2, sightT[ww2, 2])] <- 1

    m0         <- getMProbMat(gender[www], state0, rnow, monNow, ja)  #Pr t-1 to t
    wsite      <- which(sightT[, 2] > 0)
    m0[wsite, ] <- 0
    m0[cbind(wsite, sightT[wsite, 2])] <- 1

    month <- monNow
    if (month == 12) month <- 1
    m1 <- getMProbMat(gender[www], rnow, state2, month, ja)  #Pr t to t+1

    notSight <- exp(matrix(-lamda[www], nww, nreg) * matrix(effort, nww, nreg, byrow = T)) #
    wns      <- which(sightT[, 2] > 0)
    notSight[wns, ] <- 0
    notSight[cbind(wns, sightT[wns, 2])] <- 1

    h0  <- health[www, t - 1]
    wna <- which(is.na(h0))
    if (length(wna) > 0) h0[wna] <- health[www[wna], t]  #if no health at t-1, use t


    slopeSurvNow0    <- matrix(h0, nww, nreg) * bg[1]              #t-1, t
    slopeSurvNow1    <- matrix(health[www, t], nww, nreg) * bg[1]  #t, t+1
    slopeSurvNew0    <- matrix(h0, nww, nreg) * bnew[1]            #with proposed coefficients
    interceptSurvNow <- matrix(bg[-1], nww, nreg, byrow = T)
    interceptSurvNew <- matrix(bnew[-1], nww, nreg, byrow = T)

    surv0    <- inv.logit(slopeSurvNow0 + interceptSurvNow) #regional intercepts
    surv1    <- inv.logit(slopeSurvNow1 + interceptSurvNow)
    surv0[is.na(surv0)] <- surv1[is.na(surv0)]


    tmp <- moveSurv(ws0, ws1, ws2, ws3, surv0, surv1, state0, state1, state2, m0, m1, notSight)
    wsurv <- tmp$surv
    wstat <- tmp$states # check that this is the current region
    if(nww != length(wstat)){
	    tmp <- moveSurv(ws0, ws1, ws2, ws3, surv0, surv1, state0, state1, state2, m0, m1, notSight)
    	wsurv <- tmp$surv
    	wstat <- tmp$states
    }

    if (FixSurv) {
      wsurv[wsurv == 0] <- 1
    }

    deadwww <- death[www]
    deadwww[wsurv == 0] <- t
    deadwww[wsurv == 1 & deadwww == t] <- t + 1
    death[www] <- deadwww

    if (t < nt) {
      rs[www[death[www] == t], (t + 1):nt, ] <- 0
    }

    if(nww != length(wstat))stop("error in updateStates")

    rnew <- rt1*0
    rnew[cbind(c(1:nww), wstat)] <- 1
    rs[www, t, ] <- rnew

    if (t > 1) {
      wherelast <- rowMaxs(rt0 * iByr[www, ])
      wherenow  <- rowMaxs(rnew * iByr[www, ])
      tmp <- moveSum(wherelast, wherenow, gender[www], ja)
      moveFem[, , monNow]  <- moveFem[, , monNow] + tmp$fem
      moveMale[, , monNow] <- moveMale[, , monNow] + tmp$male
      moveX[, , monNow]    <- moveX[, , monNow] + tmp$x
      moveJ[, , monNow]    <- moveJ[, , monNow] + tmp$juv
    }

    survNow <- inv.logit(slopeSurvNow0 + interceptSurvNow)
    survNew <- inv.logit(slopeSurvNew0 + interceptSurvNew)

    aindex <- cbind(c(1:nww), wstat)
    dindex <- aindex[wsurv == 0, ]
    aindex <- aindex[wsurv == 1, ]

    if(length(aindex) > 0){
      psnow <- psnow + sum(log(survNow[aindex]),na.rm=T)
      psnew <- psnew + sum(log(survNew[aindex]),na.rm=T)
    }
    if(length(dindex) > 0){
      psnow <- psnow + sum(log(1 - survNow[dindex]),na.rm=T)
      psnew <- psnew + sum(log(1 - survNew[dindex]),na.rm=T)
    }
  }

if(multExp){
  movePost <- getMovePostME(maleMove, moveMale, femMove, moveFem, xMove, moveX, juvMove, moveJ)
} else {
  movePost <- getMovePost(maleMove, moveMale, femMove, moveFem, xMove, moveX, juvMove, moveJ) 
}
 


  psnow <- psnow + mydmvnorm(t(bg), priorBSurv, priorVS, log = T)
  psnew <- psnew + mydmvnorm(t(bnew), priorBSurv, priorVS, log = T)
  a <- exp(psnew - psnow)
  z <- runif(1, 0, 1)
  if (z < a) bg <- bnew

  list(rs = rs,
       death = death,
       mm = movePost$maleMoveProb,
       mf = movePost$femMoveProb,
       mx = movePost$xMoveProb,
       mj = movePost$juvMoveProb,
       bg = bg)
}
updateStates <- cmpfun(updateStates)

getMovePost <- function(maleMove, moveMale, femMove, moveFem, xMove, moveX, juvMove, moveJ){

  # This function returns movement posteriors given priors and currently imputed data, e.g.
  # maleMove is the prior
  # moveMale is the data

  for(k in 1:12){
    di <- matrix(rgamma(nreg * nreg, shape = maleMove[, , k] +
                          moveMale[, , k], scale = 1), nreg, nreg)
    maleMoveProb[, , k] <- di / matrix(colSums(di), nreg, nreg, byrow = T)

    di <- matrix(rgamma(nreg * nreg, shape = femMove[, , k] +
                          moveFem[, , k], scale = 1), nreg, nreg)
    femMoveProb[, , k] <- di / matrix(colSums(di), nreg, nreg, byrow = T)

    di <- matrix(rgamma(nreg * nreg, shape = xMove[, , k] +
                          moveX[, , k], scale = 1), nreg, nreg)
    xMoveProb[, , k] <- di / matrix(colSums(di), nreg, nreg, byrow = T)

    di <- matrix(rgamma(nreg * nreg, shape = juvMove[, , k] +
                          moveJ[, , k], scale = 1), nreg, nreg)
    juvMoveProb[, , k] <- di / matrix(colSums(di), nreg, nreg, byrow = T)
  }
  list(maleMoveProb = maleMoveProb, femMoveProb = femMoveProb, 
       xMoveProb = xMoveProb, juvMoveProb = juvMoveProb)
}

getMovePostME <- function(maleMove, moveMale, femMove, moveFem, xMove, moveX, juvMove, moveJ){
  
  # This function returns movement posteriors given priors and currently imputed data:
  # maleMove is the prior
  # moveMale is the data
  # this is for when we use multiple experts
  library(eliciteg)
  
  for(k in 1:12){
    if(k == 12) { # because right now we only have priors for Dec-->Jan moves
      
      ############ Calculate C and K for Males & Females ######################
      # Females
      cmatOutF <- matrix(NA, nrow = length(priorList$females), 
                        ncol = dim(moveFem[, , k])[2],
                        dimnames = list(c('exp1', 'exp2', 'exp3', 'exp4', 
                                          'exp5', 'exp6', 'exp7', 'exp8'),
                                        colnames(moveFem)))
      
      for(i in seq_along(1:dim(moveFem[, , k])[2])){
        for(j in seq_along(1:length(priorList$females))){
          priorTF <- priorList$females[[j]]
          cmatOutF[j, i] <- calcC(moveFem[, i, k], priorTF[, i])
        }
      }
      
      cmatOutFexp <- transformC(cmatOutF)
      kmatF <- calcK(cmatOutFexp)
      
      # Males
      cmatOutM <- matrix(NA, nrow = length(priorList$males), 
                         ncol = dim(moveMale[, , k])[2],
                         dimnames = list(c('exp1', 'exp2', 'exp3', 'exp4', 
                                           'exp5', 'exp6', 'exp7', 'exp8'),
                                         colnames(moveMale)))
      
      for(i in seq_along(1:dim(moveMale[, , k])[2])){
        for(j in seq_along(1:length(priorList$males))){
          priorTM <- priorList$males[[j]]
          cmatOutM[j, i] <- calcC(moveMale[, i, k], priorTM[, i])
        }
      }
      
      cmatOutMexp <- transformC(cmatOutM)
      kmatM <- calcK(cmatOutMexp)
      ############ End Calculate C and K for Males & Females ##################
      
      # Old Sampling -- Females
      # di <- matrix(rgamma(nreg * nreg, shape = femMove[, , k] +
      #                       moveFem[, , k], scale = 1), nreg, nreg)
      # femMoveProb[, , k] <- di / matrix(colSums(di), nreg, nreg, byrow = T)  
      
      # New Sampling -- Females
      for(j in seq_along(1:nreg)){
        
        idx <- which(rmultinom(1, 1, kmatF[, j]) == 1) # chooses the expert
        prior <- priorList$females[[idx]][, j] # gets prior for that expert
        
        di <- matrix(rgamma(nreg * 1, 
                            shape = moveFem[, j, k] + prior, 
                            scale = 1), nreg, 1) 
        
        femMoveProb[, j, k] <- di / matrix(colSums(di), nreg, 1, byrow = T) 
        
      } 
      
      # Old Sampling -- Males
      # di <- matrix(rgamma(nreg * nreg, shape = maleMove[, , k] +
      #                       moveMale[, , k], scale = 1), nreg, nreg)
      # maleMoveProb[, , k] <- di / matrix(colSums(di), nreg, nreg, byrow = T)
      
      # New Sampling -- Males
      for(j in seq_along(1:nreg)){
        
        idx <- which(rmultinom(1, 1, kmatM[, j]) == 1) # chooses the expert
        prior <- priorList$males[[idx]][, j] # gets prior for that expert
        
        di <- matrix(rgamma(nreg * 1, 
                            shape = moveMale[, j, k] + prior, 
                            scale = 1), nreg, 1) 
        
        maleMoveProb[, j, k] <- di / matrix(colSums(di), nreg, 1, byrow = T) 
        
      } 
      
    }
    
    di <- matrix(rgamma(nreg * nreg, shape = xMove[, , k] +
                          moveX[, , k], scale = 1), nreg, nreg)
    xMoveProb[, , k] <- di / matrix(colSums(di), nreg, nreg, byrow = T)
    
    di <- matrix(rgamma(nreg * nreg, shape = juvMove[, , k] +
                          moveJ[, , k], scale = 1), nreg, nreg)
    juvMoveProb[, , k] <- di / matrix(colSums(di), nreg, nreg, byrow = T)
  }
  
  list(maleMoveProb = maleMoveProb, femMoveProb = femMoveProb, 
       xMoveProb = xMoveProb, juvMoveProb = juvMoveProb)
}

getCmat <- function(cmat,br){

  c2 <- cmat[,2]
  if(c2[1] > -.5)c2[1] <- -.5
  if(c2[nrow(cmat)] > -.01)c2[nrow(cmat)] <- -.01

  dc   <- diff(c2)
  wc   <- which(dc < 0)

  while(length(wc) > 0){
    c2[wc[1]+1] <- min(c2[wc[1]] + .01,-.01)
    dc   <- diff(c2)
    wc   <- which(dc < 0)
  }

  emat <- matrix(exp(cmat[,1] + cmat[,2]*br),length(br),length(br),byrow=T)
  smat <- matrix(1,length(br),length(br))
  smat[row(smat) == col(smat)] <- 0
  c <- log(1 + rowSums(emat*smat))

  c2 <- (c - cmat[,1])/br
  if(c2[1] > -.5)c2[1] <- -.5

  c1 <- c - c2*br

  cbind(c1,c2)
}
getCmat <- cmpfun(getCmat)

proposeBreak <- function(cmat,br,brLims){

  c1 <- cmat[,1]
  c2 <- cmat[,2]

  nd <- nrow(cmat)
  brNew <- tnorm(nd,brLims[1,],brLims[2,],br,.1)

  if(nd == 1){
    c1New <- tnorm(1,1,100,cmat[1],.1)
    c2New <- (log(.5) + log(1 + exp(c1New + cmat[2]*brNew)) - c1New)/brNew
    cnew <- matrix(c(c1New,c2New),1,2)
  }

  if(nd > 1)cnew <- getCmat(cmat,brNew)

  list(cmat = cmat,br = brNew)
}
proposeBreak <- cmpfun(proposeBreak)



proposeBreakSteep <- function(breakg,steepg,bcov){

  nbreak <- length(breakg)

  minprop <- breakg - 5
  maxprop <- breakg + 5
  minprop[minprop < 5] <- 5
  maxprop[maxprop > 95] <- 95
  minprop <- c(minprop,.5)
  maxprop <- c(maxprop,nbreak)

  propvals <- tnorm.mvt(c(breakg,steepg),c(breakg,steepg),bcov,minprop,maxprop)

  propvals
}
proposeBreakSteep <- cmpfun(proposeBreakSteep)

healthPars <- function(b,priorB,priorVB,breakg,brLims,discStates,contStates,acount){

  accept <- 0

  if(acount > 10)pv1 <- log(.0005/acount)
  pv <- rlnorm(1,pv1,1)

  nbreak   <- length(breakg)

  breaks <- tnorm(nbreak,brLims[1,],brLims[2,],breakg,pv)
  db <- diff(b[,2])/2
  mids <- b[-nbreak,2] + db
  lo   <- c(b[1,2]-.2,mids)
  hi   <- c(mids,0)

  if(acount > 10)pv2 <- log(.000005/acount)
  pv <- rlnorm(1,pv2,1)

  b2 <- tnorm(nbreak,lo,hi,b[,2],pv)

  bpp <- breaks2pars(cbind(b[,1],b2),breaks)

  nss <- nbreak*2
  if(acount > 20){
    sss <- sample(1:nss,2)
    bpp[sss] <- b[sss]
  }
  if(acount > 50){
    sss <- sample(1:nss,nbreak)
    bpp[sss] <- b[sss]
  }
  if(acount > 100){
    sss <- sample(1:nss,nss-1)
    bpp[sss] <- b[sss]
  }

  discStates[is.na(discStates)] <- 0

  ww <- which(discStates > 0 & is.finite(contStates),arr.ind=T)

  pnow <- multiLogitStates(b,discStates[ww],contStates[ww],nbreak+1)
  pnew <- multiLogitStates(bpp,discStates[ww],contStates[ww],nbreak+1)

  ll1 <- sum(pnow,na.rm=T)
  ll2 <- sum(pnew,na.rm=T)

  pnow <- ll1 + mydmvnorm(as.vector(b),as.vector(priorB),priorVB,log=T)
  pnew <- ll2 + mydmvnorm(as.vector(bpp),as.vector(priorB),priorVB,log=T)

  z <- runif(1,0,1)
  r <- exp(pnew - pnow)

  like <- ll1

  if(z < r){
    b <- bpp
    breakg <- breaks
    like = ll2
    accept <- 1
  }
  list(bgg = b, breakg = breakg, accept = accept, loglik = like)
}
healthPars <- cmpfun(healthPars)

predLogit <- function(h,gchain,priorB){               #plot logit from MCMC chains

  nh   <- length(h)
  nd   <- ncol(gchain)/2
  nsim <- nrow(gchain)
  lo <- hi <- mid <- hprior <- matrix(0,nh,(nd+1))

  for(i in 1:nh){

    st <- rep(0,nsim)
    sp <- 0

    for(k in 1:nd){

      eh <- exp(gchain[,k] + gchain[,(k+nd)]*h[i])
      eh[eh > 1e+5] <- 1e+5
      tk <- eh/(1 + eh) - st
      st <- st + tk
      lo[i,k]  <- quantile(tk,.025)
      hi[i,k]  <- quantile(tk,.975)
      mid[i,k] <- quantile(tk,.5)

      ep <- exp(priorB[k,1] + priorB[k,2]*h[i])
      ep[ep > 1e+5] <- 1e+5
      tp <- ep/(1 + ep) - sp
      sp <- sp + tp
      hprior[i,k] <- tp

    }
    lo[i,nd+1]  <- 1 - quantile(st,.025)
    hi[i,nd+1]  <- 1 - quantile(st,.975)
    mid[i,nd+1] <- 1 - quantile(st,.5)

    hprior[i,nd+1] <- 1 - sp
  }

  plot(h, mid[, k], type = 'l', lwd = 1, ylim = c(0, 1), las = 1, ylab = 'Probability', xlab = 'True health')

  library(RColorBrewer)
  my.pal <- brewer.pal(3, 'Set1')
  for(k in 1:(nd+1)){

    lines(h, mid[, k], lwd = 1, col = my.pal[k])
    lines(h, lo[, k], lty = 2, col = my.pal[k])
    lines(h, hi[, k], lty = 2, col = my.pal[k])

    lines(h, hprior[, k], lty = 3, , lwd = 1, col = my.pal[k])

  }
  abline(h = .5, lty = 2)

  list(lo = lo, hi = hi, mid = mid)
}

predLogitPlot <- function(h,gchain,priorB){               #plot logit from MCMC chains with thicker lines

  nh   <- length(h)
  nd   <- ncol(gchain)/2
  nsim <- nrow(gchain)
  lo <- hi <- mid <- hprior <- matrix(0,nh,(nd+1))

  for(i in 1:nh){

    st <- rep(0,nsim)
    sp <- 0

    for(k in 1:nd){

      eh <- exp(gchain[,k] + gchain[,(k+nd)]*h[i])
      eh[eh > 1e+5] <- 1e+5
      tk <- eh/(1 + eh) - st
      st <- st + tk
      lo[i,k]  <- quantile(tk,.025)
      hi[i,k]  <- quantile(tk,.975)
      mid[i,k] <- quantile(tk,.5)

      ep <- exp(priorB[k,1] + priorB[k,2]*h[i])
      ep[ep > 1e+5] <- 1e+5
      tp <- ep/(1 + ep) - sp
      sp <- sp + tp
      hprior[i,k] <- tp

    }
    lo[i,nd+1]  <- 1 - quantile(st,.025)
    hi[i,nd+1]  <- 1 - quantile(st,.975)
    mid[i,nd+1] <- 1 - quantile(st,.5)

    hprior[i,nd+1] <- 1 - sp
  }

  plot(h, mid[, k], type = 'l', lwd = 2.5, ylim = c(0, 1), las = 1, ylab = 'Probability', xlab = 'True health')

  library(RColorBrewer)
  my.pal <- brewer.pal(3, 'Set1')
  for(k in 1:(nd+1)){

    lines(h, mid[, k], lwd = 2.5, col = my.pal[k])
    lines(h, lo[, k], lty = 2, lwd = 1.5,col = my.pal[k])
    lines(h, hi[, k], lty = 2, lwd = 1.5,col = my.pal[k])

    lines(h, hprior[, k], lty = 3, , lwd = 1, col = my.pal[k])

  }
  abline(h = .5, lty = 2)

  list(lo = lo, hi = hi, mid = mid)
}


plotLogit <- function(lims,breaks,h,cmat){

  tmp  <- pars2p(cmat,h)
  library(RColorBrewer)
  my.pal <- brewer.pal(ncol(lims) + 1, 'Set1')

  plot(h, tmp[, 1], type = 'l', lwd = 2, xlab = 'Latent health scale', ylab = 'Probabilty', ylim = c(0, 1))

  for(j in 1:ncol(lims)){
    polygon(c(lims[, j], rev(lims[, j])), c(0, 0, 1, 1), border = my.pal[j])
    lines(h, tmp[, j], col = my.pal[j], lwd = 2)
    l1   <- 0
    if(j > 1)l1 <- breaks[j - 1]
    midx <- (l1 + breaks[j])/2
    #     text(midx, 0.8, j,col = my.pal[j], cex = 1.4)
  }
  lines(h, tmp[, j + 1], col = my.pal[j + 1], lwd = 2)
  midx <- ( breaks[j] + 100)/2
  #   text(midx, 0.8, j + 1, col = my.pal[j + 1], cex = 1.4)
  abline(h = 0.5, lty = 2)
}


healthStateSpacePars <- function(){   #these are bq

  x1 <- firstSight
  x2 <- death - 1
  y1 <- firstSight + 1
  y2 <- death

  nx <- sum(x2 - x1) + n
  xx <- matrix(1,nx,5)
  yy <- rep(0,nx)

  for(i in 1:n){

    xi <- x1[i]:x2[i]
    yi <- y1[i]:y2[i]

    if(max(yi) > nt){
      wi <- which(yi <= nt)
      xi <- xi[wi]
      yi <- yi[wi]
    }

    age <- xi/12
    ja  <- age*0
    ja[age <= juvAge] <- 1

    xx[xi, 2] <- health[i, xi]
    xx[xi, 3] <- ja
    xx[xi, 4] <- age
    xx[xi, 5] <- age^2
    yy[yi] <- health[i, yi]

  }

  #remove individuals alive once:
  wna <- which(is.na(xx[,2]) | is.na(yy))
  if(length(wna) > 0){
    xx <- xx[-wna,]
    yy <- yy[-wna]
  }

  V <- solve(crossprod(xx)/herror + priorIVQ)
  v <- crossprod(xx,yy)/herror + priorIVQ %*% priorBQ
  mu <- V%*%v

  bq <- tnorm.mvt(bq, mu, V, loBQ, hiBQ)
  if(!JUV){bq[varQ == 'juv'] <- 0}
  s1 <- h1 + .5 * nx
  s2 <- h2 + .5 * sum( (yy - xx %*% bq) ^ 2 )
  herror <- 1/rgamma(1, s1, s2)
  list(bq = bq, herror = herror)
}
healthStateSpacePars <- cmpfun(healthStateSpacePars)


healthStateSpacePars2 <- function(){   #these are bq; fixes possible indexing error for xi

  x1 <- firstSight
  x2 <- death - 1
  y1 <- firstSight + 1
  y2 <- death

  nx <- sum(x2 - x1) + n
  xx <- matrix(1,nx,5)
  yy <- rep(0,nx)
  counter <- 0
  for(i in 1:n){

    xi <- x1[i]:x2[i]
    yi <- y1[i]:y2[i]

    if(max(yi) > nt){
      wi <- which(yi <= nt)
      xi <- xi[wi]
      yi <- yi[wi]
    }

#     age <- xi/12
    age <-  (xi - firstSight[i])/12 + 1
    ja  <- age*0
    ja[age <= juvAge] <- 1

    xx[counter + (1:length(xi)), 2] <- health[i, xi]
    xx[counter + (1:length(xi)), 3] <- ja
    xx[counter + (1:length(xi)), 4] <- age
    xx[counter + (1:length(xi)), 5] <- age^2
    yy[counter + (1:length(yi))] <- health[i, yi]

    counter <- counter + length(xi)
  }
  xx = xx[1:counter,]
  yy = yy[1:counter]

  #remove individuals alive once:
  wna <- which(is.na(xx[,2]) | is.na(yy))
  if(length(wna) > 0){
    xx <- xx[-wna,]
    yy <- yy[-wna]
  }

  V <- solve(crossprod(xx)/herror + priorIVQ)
  v <- crossprod(xx,yy)/herror + priorIVQ %*% priorBQ
  mu <- V%*%v


  bq <- tnorm.mvt(bq, mu, V, loBQ, hiBQ)
  if(!JUV){bq[varQ == 'juv'] <- 0}
  s1 <- h1 + .5 * nx
  s2 <- h2 + .5 * sum( (yy - xx %*% bq) ^ 2 )
  herror <- 1/rgamma(1, s1, s2)
  list(bq = bq, herror = herror)
}
healthStateSpacePars2 <- cmpfun(healthStateSpacePars2)

multiLogitStates <- function(b, discStates, contStates, maxD){    # update of continuous states based on discrete states

  tiny <- 1e-20
  nn   <- length(discStates)

  discStates[is.na(discStates)] <- 0
  prob <- discStates*0 + 1

  if(maxD == 1)return(log(prob))

  wk <- which(discStates == 1 & is.finite(contStates), arr.ind = T)

  prob[wk] <- invlogt(b[1, ], contStates[wk])

  if(maxD > 2){
    for(k in 2:(maxD-1)){
      wk       <- which(discStates == k & is.finite(contStates), arr.ind = T)
      prob[wk] <- invlogt(b[k, ], contStates[wk]) -
        invlogt(b[(k - 1), ], contStates[wk])
    }
  }
  prob[discStates == maxD] <- 1 - invlogt(b[(maxD - 1), ], contStates[discStates == maxD])

  prob[prob < tiny] <- tiny
  log(prob)
}


propH <- function(cnow,hoffset = 5){

  clo <- cnow - hoffset
  clo[clo < 0] <- 0
  chi <- cnow + hoffset
  chi[chi > 100] <- 100

  tnorm(length(cnow), clo, chi, cnow, rexp(length(cnow), 1/10))

}

multiLogitStates <- cmpfun(multiLogitStates)



updateHealth <- function(){

  bh <- matrix(bh,2,2)
#   bSkin <- matrix(bSkin,2,2)

  accept <- 0
  total  <- 0
  hsNew  <- hStates*0
  skinNew<- skin*0

  health[health < .1]  <- .1
  hStates[hStates < 1] <- 1

  # hoffset <- 5

  t1 <- min(firstSight)

  for(t in t1:nt){

    health[death < t,t] <- 0
    wa <- which(firstSight <= t & death >= t)
    na <- length(wa)
    if(na == 0)next

    age <-  (t - firstSight[wa])/12 + 1
    ja  <- age * 0
    ja[age <= juvAge] <- 1

    hnow  <- hStates[wa, t]    #current body fat states

    hlast <- hnow
    hnext <- hnow

    if(t > 1) hlast <- hStates[wa, t - 1]
    if(t < nt)hnext <- hStates[wa, t + 1]

    hnow[hnow == 0] <- hlast[hlast == 0] <- hnext[hnext == 0] <- 1

    cnow <- health[wa, t]
    cnew <- propH(cnow, hoffset)

    if(bfFlag){pnow <- multiLogitStates(bh, hStates[wa, t], cnow, 3)
               pnew <- multiLogitStates(bh, hStates[wa, t], cnew, 3)}

    if(entFlag){pnow <- pnow + multiLogitStates(bEnt, tStates[wa, t], cnow, length(breakEnt) + 1)
                pnew <- pnew + multiLogitStates(bEnt, tStates[wa, t], cnew, length(breakEnt) + 1)}

    if(pFlag){pnow <- pnow + multiLogitStates(bProp, pStates[wa, t], cnow, length(breakProp) + 1)
              pnew <- pnew + multiLogitStates(bProp, pStates[wa, t], cnew, length(breakProp) + 1)}

    if(rkFlag){pnow <- pnow + multiLogitStates(bRake, rake[wa, t], cnow, length(breakRake) + 1)
               pnew <- pnew + multiLogitStates(bRake, rake[wa, t], cnew, length(breakRake) + 1)}

    if(skFlag){pnow <- pnow + multiLogitStates(bSkin, skin[wa, t], cnow, length(breakSkin) + 1)
               pnew <- pnew + multiLogitStates(bSkin, skin[wa, t], cnew, length(breakSkin) + 1)}

    if(cyFlag){pnow <- pnow + multiLogitStates(bCyam, cyam[wa, t], cnow, length(breakCyam) + 1)
               pnew <- pnew + multiLogitStates(bCyam, cyam[wa, t], cnew, length(breakCyam) + 1)}

    if(caFlag){pnow <- pnow + multiLogitStates(bCalf, calves[wa, t], cnow, length(breakCalf) + 1)
               pnew <- pnew + multiLogitStates(bCalf, calves[wa, t], cnew, length(breakCalf) + 1)}

    hnow  <- health[wa, t]
    hlast <- hnow
    hnext <- hnow
    if(t > t1)hlast <- health[wa, t - 1]
    if(t < nt)hnext <- health[wa, t + 1]

    #survival

    if(t < nt){

      w2 <- which(death[wa] > (t+1) & is.finite(hlast))  #does not die before t+1
      wd <- which(death[wa] %in% t:(t+1))                    #dies before t+1

      co <- matrix(0,na,nreg+1)
      co[,-1] <- rState[wa,t,]
      co[,1]  <- cnow
      ce      <- co
      ce[,1]  <- cnew

      xnow <- inv.logit(co%*%bg)
      xnew <- inv.logit(ce%*%bg)
      xnow[xnow > .9999999] <- .9999999
      xnew[xnew > .9999999] <- .9999999

      pnow[w2] <- pnow[w2] + log(xnow[w2])
      pnew[w2] <- pnew[w2] + log(xnew[w2])

      pnow[wd] <- pnow[wd] + log(1 - xnow[wd])
      pnew[wd] <- pnew[wd] + log(1 - xnew[wd])
    }

    #process
    V <- rep(0,na)
    v <- V
    z0 <- matrix(0,na,length(varQ))
    z1 <- z0

    w0 <- which(death[wa] > t & is.finite(hlast))      #present at t-1
    w2 <- which(death[wa] > (t+1) & is.finite(hnext))  #present at t+1

    if(length(w0) > 0){
      z0[w0, 1] <- 1
      z0[w0, 2] <- hlast[w0]
      z0[w0, 3] <- ja[w0]
      z0[w0, 4] <- age[w0]
      z0[w0, 5] <- age[w0]^2
      mu <-  z0[w0, ] %*% bq
      pnow[w0] <- pnow[w0] + dnorm(cnow[w0], mu, sqrt(herror), log = T)
      pnew[w0] <- pnew[w0] + dnorm(cnew[w0], mu, sqrt(herror), log = T)
    }
    if(length(w2) > 0){
      z1[w2,1] <- 1
      z1[w2,2] <- hnow[w2]
      z1[w2,3] <- ja[w2]
      z1[w2,4] <- age[w2] + 1/12
      z1[w2,5] <- (age[w2] + 1/12)^2
      mu1 <-  z1[w2,]%*%bq
      z1[w2,2] <- cnew[w2]
      mu2 <- z1[w2,]%*%bq
      pnow[w2] <- pnow[w2] + dnorm(hnext[w2],mu1,sqrt(herror),log=T)
      pnew[w2] <- pnew[w2] + dnorm(hnext[w2],mu2,sqrt(herror),log=T)
    }


    a <- exp(pnew - pnow)
    z <- runif(na, 0, 1)
    cnow[z < a] <- cnew[z < a]
    health[wa, t] <- cnow
    if(t < nt)health[death < t, (t + 1):nt] <- .1
    total <- total + na
    accept <- accept + length(z[z < a])

    if(missDataImp){
      # Body Fat
#       hStatePrior[wa, t][is.na(hStatePrior[wa, t])] <- 3
      hProb <- predProb(bh, cnow) * t( hTable[, hStatePrior[wa, t]] ) # old code: predProb(bh, cnow) * hStatePrior[wa, ]
      hProb[hProb < 0] <- 0
      hProb <- hProb/matrix(rowSums(hProb), na, 3)
      hh <- rowSums(myrmultinom(1, hProb) * matrix(c(1:3), na, 3, byrow = T))
      hsNew[wa, t] <- hh

      # Skin
      skProb <- predSkinProb(bSkin, cnow) * t( skTable[, skinStatePrior[wa, t]] )
      skProb[skProb < 0] <- 0
      skProb <- skProb/matrix(rowSums(skProb), na, 2)
      sk <- rowSums(myrmultinom(1, skProb) * matrix(c(1:2), na, 2, byrow = T))
      skinNew[wa, t] <- sk
    }

  }

  if(missDataImp){
    hStates[missH] <- hsNew[missH]
    skin[missSkin] <- skinNew[missSkin]
  }


  list(health = health, hStates = hStates, skin = skin, accept = accept/total)
}
updateHealth <- cmpfun(updateHealth)

myrmultinom <- function(size,p){

  #n multinomial r.v. for a n by ncol(p) matrix of probs
  #each row of p is a probability vector

  n     <- nrow(p)
  J     <- ncol(p)
  y     <- matrix(0,n,J)
  sizej <- size
  if(length(size) == 1)sizej <- rep(size,n)
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
    pj    <- matrix(p[,c((j+1):J)]/dpj,n,J-j)
    wj    <- which(sumj < size,arr.ind=T)
  }

  y[,J] <- size - rowSums(y)
  y
}
myrmultinom <- cmpfun(myrmultinom)


predProb <- function(b,h){

  pn  <- matrix(0,length(h),3)
  pn[,1] <- invlogt(b[c(1,3)],h)
  pn[,2] <- invlogt(b[c(2,4)],h) - invlogt(b[c(1,3)],h)
  pn[,3] <- 1 - invlogt(b[c(2,4)],h)
  pn
}

predSkinProb <- function(b,h){

  pn  <- matrix(0,length(h), 2)
  pn[,1] <- invlogt(b, h)
  pn[,2] <- 1 - invlogt(b, h)
  pn
}


invlogt <- function(pc,lss){
  z <- exp(pc[1] + pc[2]*lss)
  z/(1 + z)
}

invlogt <- cmpfun(invlogt)


priorHC3 <- function(c11,halves,steep){ #halves - where Pr = 1/2

  c01 <- -c11*halves[1]
  c12 <- c01/halves[2] - steep
  c02 <- -c12*halves[2]

  hseq <- c(0:100)
  cmat <- matrix(c(c01,c02,c11,c12),2,2)  #c11 < 0
  hout <- predH(cmat,hseq)

  plot(hout[,1],type='l',ylim=c(0,1))
  lines(hout[,2],col=2)
  lines(hout[,3],col=3)

  abline(h=.5,lty=2)
  abline(v=halves,lty=10)

  cmat
}
priorHC3 <- cmpfun(priorHC3)

pars2p <- function(cmat,h){

  nd <- nrow(cmat)
  nh <- length(h)

  c1 <- matrix(cmat[,1],nh,nd,byrow=T)
  c2 <- matrix(cmat[,2],nh,nd,byrow=T)
  hh <- matrix(h,nh,nd)
  eh <- exp(c1 + c2*hh)

  theta <- matrix(0,nh,nd+1)
  sumt  <- rep(0,nh)

  for(k in 1:nd){

    tk   <- eh[,k]/(1 + eh[,k]) - sumt
    theta[,k] <- tk
    sumt <- sumt + tk
  }
  theta[,nd+1] <- 1 - rowSums(theta)

  theta
}
pars2p <- cmpfun(pars2p)


makePriorB <- function(ir, sr, breaks, varInt = 10, varSlp = 10){

  n  <- length(breaks)
  mm <- matrix(cbind(seq(ir[1], ir[2], length.out = n), seq(sr[1], sr[2], length.out = n)), n, 2)
  cmat <- breaks2pars(mm, breaks)
  tmp  <- pars2p(cmat, hseq)

  pvar <- diag(c(rep(varInt, n), rep(varSlp, n)))

  list(meanPrior = cmat, varPrior = pvar)
}
makePriorB <- cmpfun(makePriorB)

breaks2pars <- function(cmat,breaks){

  nd <- length(breaks)
  c0 <- rep(0,nrow(cmat))

  for(k in 1:nd){

    qk <- 0
    if(k > 1){
      for(j in 1:(k-1)){
        ej <- exp(c0[j] + cmat[j,2]*breaks[k])
        qk <- qk + ej/(1 + ej)
      }
    }
    D     <- -log(1/(.5 + qk) - 1)
    cc0   <- D - cmat[k,2]*breaks[k]
    if(k > 1){
      if(cc0 < cmat[k-1,1]){
        cc0 <- cmat[k-1,1] + .5
        cmat[k,2] <- (D - cc0)/breaks[k]
      }
    }
    cmat[k,1] <- cc0
  }
  cmat
}
breaks2pars <- cmpfun(breaks2pars)

breaks2mids <- function(breaks,sss){

  bb <- c(0,breaks) + diff(c(0,breaks,100))/2
  bb[sss]

}

breaks2mids <- cmpfun(breaks2mids)



priorHC <- function(c11,halves,steep){ #halves - where Pr = 1/2

  nd <- length(halves)
  cc <- matrix(0,nd,2)
  cc[1,2] <- c11

  for(k in 1:(nd-1)){
    cc[k,1]   <- -cc[k,2]*halves[k]

    if(k > 1){
      if(cc[k,1] <= cc[k-1,1]) cc[k,1] <- cc[k-1,1] + 1
    }
    cc[k+1,2] <- cc[k,1]/halves[k+1] - steep
    if(cc[k+1,2] <= cc[k,2])cc[k+1,2] <- cc[k,2] + .5/length(halves)
    cc[k+1,1] <- -cc[k+1,2]*halves[k+1]
    if(cc[k+1,1] <= cc[k,1])cc[k+1,1] <- cc[k,1] + 1
  }

  if(cc[nd,1] == cc[nd-1,1])cc[nd,1] <- cc[nd,1] + .5/length(halves)

  cc
}
priorHC <- cmpfun(priorHC)

InitialMove <- function(df,wt){
  # df will be the narrative priors as a data frame
  #next line just for testing
  # df <- adultMales

  move <- array(0,dim=c(nreg,nreg,12))
  dimnames(move)[[3]] <- c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
  dimnames(move)[[1]] <- dimnames(move)[[2]] <- regID
  moveProb <- move

  monthvec <- c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
  for(k in 1:12){
    df.sub <- df[df$month == monthvec[k],]
    df.sub <- df.sub[,c(4:13)]
    rownames(df.sub) <- colnames(df.sub)

    move[,,k] <- as.matrix(df.sub[match(regID, rownames(df.sub)), match(regID, colnames(df.sub))])
    moveProb[,,k] <- move[,,k] / matrix(colSums(move[,,k]), nreg,nreg,byrow=T)
  }

  maxmove <- max(move)
  minmove <- min(move)

  move    <- wt*(move - minmove)/(maxmove - minmove) + 1

  list(moveProb = moveProb, move = move)
}


InitialMoveOLD<- function(){

  moveProb <- nextmat[match(regID,rownames(nextmat)),match(regID,colnames(nextmat))]

  moveProb[moveProb == 0] <- .0001
  msum <- apply(moveProb,2,sum)
  moveProb <- moveProb/matrix(msum,nreg,nreg,byrow=T)
  mWt <- 100
  bWt <- 10

  #jan movement
  janm <- moveProb*0 + bWt
  janm[c('GOM','NE'),'GOM']  <- mWt
  janm[c('NE','GSC','BOF','SEUS','MIDA'),]  <- 0

  janf <- janm*0 + bWt
  janf[c('SEUS','NE'),'MIDA'] <- mWt
  janf[c('MIDA','SEUS','NE'),'SEUS'] <- mWt
  janf[c('NE','GSC','BOF','RB','EAST'),] <- 0

  #feb movement
  febm <- moveProb*0 + bWt
  febm[c('GOM','NE'),'GOM']  <- mWt
  febm[c('BOF','SEUS','MIDA'),]  <- 0

  febf <- febm*0 + 1
  febf[c('NE'),'MIDA'] <- mWt
  febf[c('MIDA','SEUS','NE'),'SEUS'] <- mWt
  febf[c('BOF','RB','EAST'),] <- 0

  #mar movement
  marm <- moveProb*0 + bWt
  marm[c('GOM','NE','GSC'),'NE']  <- mWt
  marm[c('BOF','SEUS','MIDA'),]  <- 0

  marf <- marm*0 + bWt
  marf[c('GSC','NE'),'MIDA'] <- mWt
  marf[c('MIDA','GSC','GOM','NE'),'SEUS'] <- mWt
  marf[c('RB','BOF','EAST','NRTH'),] <- 0

  #apr movement
  aprm <- moveProb*0 + bWt
  aprm[c('GOM','NE'),'GOM']  <- mWt
  aprm[c('BOF','SEUS','MIDA'),]  <- 0

  aprf <- aprm*0 + bWt
  aprf[c('SEUS','NE'),'MIDA'] <- mWt
  aprf[c('MIDA','SEUS','NE'),'SEUS'] <- mWt
  aprf[c('RB','BOF','EAST','NRTH'),] <- 0

  #may movement
  maym <- moveProb*0 + bWt
  maym[c('GOM','GSC'),'NE']  <- mWt
  maym[c('SEUS','MIDA'),]  <- 0

  mayf <- maym*0 + bWt
  mayf[c('GOM','GSC'),'NE'] <- mWt
  mayf[c('MIDA','GSC','NE'),'SEUS'] <- mWt
  mayf[c('SEUS'),]  <- 0

  #jun movement
  junm <- moveProb*0 + bWt
  junm[c('GOM','GSC','BOF','JL'),'GOM']  <- mWt
  junm[c('GOM','GSC','BOF','JL'),'GSC']  <- mWt
  junm[c('SEUS','MIDA'),]  <- 0

  junf <- junm

  #jul movement
  julm <- moveProb*0 + bWt
  julm[c('BOF','RB'),c('GOM','BOF','JL')]  <- mWt
  julm[c('SEUS','MIDA'),]  <- 0

  julf <- julm*0 + bWt
  julf[c('BOF'),c('GOM','BOF','JL')] <- mWt
  julf[c('SEUS','MIDA'),]  <- 0

  #aug movement
  augm <- moveProb*0 + bWt
  augm[c('GOM','BOF','RB'),'BOF']  <- mWt
  augm[c('SEUS','MIDA'),]  <- 0

  augf <- augm*0 + bWt
  augf[c('BOF'),'BOF'] <- mWt
  augf[c('SEUS','MIDA'),]  <- 0

  #sep movement
  sepm <- moveProb*0 + bWt
  sepm[c('GOM','BOF','RB','JL'),'BOF']  <- mWt
  sepm[c('SEUS','MIDA'),]  <- 0

  sepf <- sepm*0 + bWt
  sepf[c('GOM','BOF','RB','JL','SEUS','MIDA'),'BOF'] <- mWt
  sepf[c('SEUS','MIDA'),]  <- 0

  #oct movement
  octm <- moveProb*0 + bWt
  octm[c('GOM','JL'),'GOM']  <- mWt
  octm[c('GOM','BOF','RB','JL'),'BOF']  <- mWt
  octm[c('SEUS','MIDA'),]  <- 0

  octf <- octm*0 + bWt
  octf[c('GOM','BOF','RB','JL','MIDA','SEUS'),'BOF'] <- mWt
  octf[c('SEUS','MIDA'),]  <- 0

  #nov movement
  novm <- moveProb*0 + bWt
  novm[c('GOM','JL'),'GOM']  <- mWt
  novm[c('SEUS','MIDA','BOF','GSC'),]  <- 0

  novf <- novm*0 + bWt
  novf[c('SEUS','MIDA'),'MIDA'] <- mWt
  novf[c('JL','GOM'),'GOM'] <- mWt

  #dec movement
  decm <- moveProb*0 + bWt
  decm[c('GOM'),'GOM']  <- mWt
  decm[c('SEUS','MIDA','BOF','RB','GSC'),]  <- 0

  decf <- decm*0 + bWt
  decf[c('SEUS'),'MIDA'] <- mWt
  decf[c('SEUS'),'SEUS'] <- mWt
  decf[c('EAST','BOF','RB','GSC'),]  <- 0

  maleMove <- array(0,dim=c(nreg,nreg,12))
  dimnames(maleMove)[[3]] <- c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
  dimnames(maleMove)[[1]] <- dimnames(maleMove)[[2]] <- regID

  femMove  <- maleMove
  maleMove[,,1] <- janm
  femMove [,,1] <- janf
  maleMove[,,2] <- febm
  femMove [,,2] <- febf
  maleMove[,,3] <- marm
  femMove [,,3] <- marf
  maleMove[,,4] <- aprm
  femMove [,,4] <- aprf
  maleMove[,,5] <- maym
  femMove [,,5] <- mayf
  maleMove[,,6] <- junm
  femMove [,,6] <- junf
  maleMove[,,7] <- julm
  femMove [,,7] <- julf
  maleMove[,,8] <- augm
  femMove [,,8] <- augf
  maleMove[,,9] <- sepm
  femMove [,,9] <- sepf
  maleMove[,,10] <- octm
  femMove [,,10] <- octf
  maleMove[,,11] <- novm
  femMove [,,11] <- novf
  maleMove[,,12] <- decm
  femMove [,,12] <- decf

  maleMoveProb <- femMoveProb <- xMoveProb <- femMove*0

  maleMove[maleMove == 0] <- .01
  femMove[femMove == 0] <- .01

  for(k in 1:12){

    maleMoveProb[,,k] <- maleMove[,,k]/matrix(apply(maleMove[,,k],2,sum),nreg,nreg,byrow=T)
    femMoveProb[,,k]  <- femMove[,,k]/matrix(apply(femMove[,,k],2,sum),nreg,nreg,byrow=T)
    xMoveProb[,,k] <- (femMoveProb[,,k] + maleMoveProb[,,k])/2
  }

  list(maleMoveProb = maleMoveProb, femMoveProb = femMoveProb, xMoveProb = xMoveProb,
       maleMove = maleMove, femMove = femMove, xMove = (maleMove + femMove)/2)
}


plotOne <- function(egno,sights,calves,batch){
  require(ggplot2)
  require(RColorBrewer)
  require(gridExtra)
  require(stringr)
  require(lubridate)
  dsub <- sights[which(sights$SightingEGNo==egno),]
  bsub <- batch[which(batch$EGNo==egno),]
  tsub  <- dsub[dsub$scar != 'FALSE',]

  if(nrow(tsub)>0){  events <- unique(tsub$EntglEventNo)
                     trect <- data.frame(event=events, start=ymd('1000/01/01'), end=ymd('1000/01/01'), severity=NA, gear=NA, comb=NA, scar = NA)
                     for(i in 1:length(events)){
                       trect[trect$event==events[i][1],'start'] <- as.Date(tsub[tsub$EntglEventNo==events[i],'EntglStartDate'][1])
                       trect[trect$event==events[i][1],'end'] <- as.Date(tsub[tsub$EntglEventNo==events[i],'EntglEndDate'][1])
                       trect[trect$event==events[i][1],'severity'] <- tsub[tsub$EntglEventNo==events[i],'severity'][1]
                       trect[trect$event==events[i][1],'gear'] <- tsub[tsub$EntglEventNo==events[i],'gear'][1]
                       trect[trect$event==events[i][1],'scar'] <- tsub[tsub$EntglEventNo==events[i],'scar'][1]
                     }
                     trect[is.na(trect$start),'start'] <- as.Date(min(dsub$Date))
                     trect$comb <- paste(str_sub(trect$severity,1,3),trect$gear,sep="")

                     if('POSSIBLE' %in% trect$scar){
                       ddsub <- subset(trect, scar=='POSSIBLE')
                       startline <- as.Date(ddsub$start)
                     }
  }


  cyr <- as.numeric(calves[calves$EGNo==egno,'Calving.Year'])
  if(length(cyr)>0){
    cyrf <- cyr +1
    cyrp <- cyr -1
    call <- sort(c(cyr,cyrp,cyrf))
    crect <- data.frame(year=call,interval=NA,start=ymd('1000/01/01'),end=ymd('1000/01/01'))
    crect[crect$year %in% cyr,'interval'] <- 'Calving'
    crect[crect$year %in% cyrp,'interval'] <- 'Gestation'
    crect[crect$year %in% cyrf,'interval'] <- 'Recovery'
    crect[,'start'] <- as.Date(paste(crect[,'year'],1,1,sep="/"))
    crect[,'end'] <- as.Date(paste(crect[,'year'],12,31,sep="/"))}

  gender <- 'M'
  if('F' %in% dsub$GenderCode){gender <- 'F'}
  base_size=12
  ind <-   qplot(x=as.Date(Date),0,data=dsub,geom='line') +
    geom_line(colour=grey(.75))+
    scale_x_date(major="5 years",minor="1 years",limits=as.Date(c("1990-01-01","2012-01-01")))+
    ylab("")+xlab("")+
    scale_y_discrete(breaks = NA) +
    theme_set(theme_bw())+
    labs(colour="The Data Range")+
    labs(size="The Data Range")+
    opts(axis.text.x = theme_text(size = base_size,vjust=1),
         axis.text.y = theme_text(size = base_size,hjust=1),
         axis.title.x=theme_text(size=base_size),
         axis.title.y=theme_text(size=base_size,angle=90),
         legend.title=theme_text(size=base_size),
         legend.text=theme_text(size=base_size-2))


  if(gender == 'F' & length(cyr > 0)){		ind.ts <-   ind +
                                          geom_point(data=bsub,aes(as.Date(Date),0,size=bsub$PhotoCount))+
                                          geom_rect(aes(NULL, NULL, xmin = as.Date(start), xmax = as.Date(end), fill = interval), ymin = -Inf, ymax = Inf, data = crect)+
                                          scale_fill_manual(values=alpha(c("#E41A1C", "#377EB8", "#4DAF4A"), 0.2), breaks = c('Gestation', 'Calving', 'Recovery'),
                                                            labels = c('Gestation', 'Calving', 'Recovery'))+
                                          opts(title=paste("EGNo:",egno," (",gender,')',", # of Sightings per Batch",sep=""))} else {	ind.ts <-   ind +
                                                                                                                                        geom_point(data=bsub,aes(as.Date(Date),0,size=bsub$PhotoCount))+
                                                                                                                                        opts(title=paste("EGNo:",egno," (",gender,')',", # of Sightings per Batch",sep=""))}

  if(nrow(tsub)>0){		ind.en <-   ind +
                       geom_point(data=bsub,aes(as.Date(Date),0),size=3)+
                       geom_rect(aes(NULL, NULL, xmin = as.Date(start), xmax = as.Date(end), fill = comb), ymin = -Inf, ymax = Inf, data = trect)+
                       scale_fill_manual('The Data Range',values=alpha(c('min0' =  "#A6CEE3", 'min1' = "#1F78B4",
                                                                         'mod0' =  "#B2DF8A", 'mod1' = "#33A02C", 'sev0' = "#FB9A99", 'sev1' = "#E31A1C"),0.5))+
                       opts(title="Entanglement Events")}

  if(exists('startline')){
    ind.en <- ind.en + geom_vline(aes(xintercept = startline))
  }

  ind.fat <-   ind +
    geom_point(pch=1,colour='grey50',size=3)+
    geom_point(data=subset(dsub,BodyFatCode == 1 | BodyFatCode == 2 | BodyFatCode == 3),aes(colour=factor(BodyFatCode)),size=4)+
    scale_colour_brewer(palette='Dark2',breaks=c(1,2,3),labels=c("Not Thin","Thin",'Very Thin'))+
    opts(title="Body Fat Levels")

  ind.bugs <-   ind +
    geom_point(pch=1,colour='grey50',size=2)+
    geom_point(data=subset(dsub,CyamidCode != NA),aes(colour=factor(CyamidCode)),pch=1,size=2)+
    scale_colour_brewer(palette='Dark2',	breaks=c(1,2,'X'),labels=c("Few","Many",'Not seen'))+
    opts(title="Cyamids on Blow Holes")

  ind.skin <-   ind +
    geom_point(pch=1,colour='grey50',size=2)+
    geom_point(data=subset(dsub,SkinCode != NA),aes(colour=factor(SkinCode)),pch=1,size=2)+
    scale_colour_brewer(palette='Dark2',	breaks=c(1,2,'X'),labels=c("Healthy","Poor",'Not seen'))+
    opts(title="Overall Skin Condition")

  ind.rr <-   ind +
    geom_point(pch=1,colour='grey50',size=2)+
    geom_point(data=subset(dsub,RightRakeMarkCode != NA),aes(colour=factor(RightRakeMarkCode)),pch=1,size=2)+
    scale_colour_brewer(palette='Dark2',	breaks=c(1,2,3,'X'),labels=c("Few","Moderate",'Severe','Not Seen'))+
    opts(title="Rakemarks, Right Side")

  ind.lr <-   ind +
    geom_point(pch=1,colour='grey50',size=2)+
    geom_point(data=subset(dsub,LeftRakeMarkCode != NA),aes(colour=factor(LeftRakeMarkCode)),pch=1,size=2)+
    scale_colour_brewer(palette='Dark2',	breaks=c(1,2,3,'X'),labels=c("Few","Moderate",'Severe','Not Seen'))+
    opts(title="Rakemarks, Left Side")

  # if(nrow(trect)>0){grid.arrange(ind.ts, ind.en, ind.fat, ind.bugs, ind.skin, ind.rr, ind.lr,nrow=7)}else{
  # grid.arrange(ind.ts, ind.fat, ind.bugs, ind.skin, ind.rr, ind.lr,nrow=6)}

  if(exists('trect')){grid.arrange(ind.ts, ind.en, ind.fat, nrow = 3)}else{
    grid.arrange(ind.ts, ind.fat, nrow= 2)}

  # grid.arrange(ind.bugs, ind.skin, ind.lr, ind.rr, nrow = 4)
}


payrescale<-function(oldob,payob,newob){
  oldrange<-range(oldob,na.rm=T)
  newrange<-range(newob,na.rm=T)
  payob01<-(payob-oldrange[1])/diff(oldrange)
  return(payob01*diff(newrange)+newrange[1])
}


makeHPrior <- function(hObs, hStates, skinFlag){

  hh <- hObs*0 + 1
  tvec  <- 1:ncol(hObs)
  hTobs <- hStates
  nt    <- ncol(hObs)
  n     <- nrow(hObs)


  if(skinFlag){
    hTobs[hStates == 1] <- 1 # imputed sightings
    hTobs[hStates == 2] <- 5
  } else {
    hTobs[hStates == 1] <- 1 # imputed sightings
    hTobs[hStates == 2] <- 3
    hTobs[hStates == 3] <- 5
  }

  hp <- hTobs

  for(i in 1:n){
    for(t in 1:nt){
      if(is.finite(hTobs[i,t]))next
      nlast <- max(hh[i,1:t]*tvec[1:t],na.rm=T)
      nnext <- min(hh[i,t:nt]*tvec[t:nt],na.rm=T)
      tlast  <- hTobs[i, nlast ]
      tnext  <- hTobs[i, nnext ]
      h0     <- hTobs[i,nlast]
      if(!is.finite(nnext))nnext <- t
      if(!is.finite(tlast))tlast <- tnext
      if(!is.finite(tnext))tnext <- tlast
      if(!is.finite(nlast)){
        nlast <- t
        h0    <- tlast
      }

      #       denom  <- nnext - nlast
      slope  <- (tnext - tlast)/(nnext - nlast)
      hp[i,t] <- h0 + slope*(t - nlast)
    }
  }

  round(hp,0)
}

makeHPrior6months <- function(hObs, hStates, skinFlag){
  # Unlike makeHPrior, the goal of this is to only interpolate when we are +/- 6 months
  # of a sighting that includes values for the visual health parameters

  hh <- hObs*0 + 1
  tvec  <- 1:ncol(hObs)
  hTobs <- hStates
  nt    <- ncol(hObs)
  n     <- nrow(hObs)

  hTobs <- hStates
  if(skinFlag){
    hTobs[hStates == 1] <- 1 # imputed sightings
    hTobs[hStates == 2] <- 5
  } else {
    hTobs[hStates == 1] <- 1 # imputed sightings
    hTobs[hStates == 2] <- 3
    hTobs[hStates == 3] <- 5
  }

  hp <- hTobs

  for(i in 1:n){
    firstObs <- min(which(is.finite(hTobs[i,])), na.rm = T)
    lastObs <- max(which(is.finite(hTobs[i,])), na.rm = T)

    for(t in 1:nt){
      if(is.finite(hTobs[i, t]))next

#       if(t < firstObs & abs(firstObs - t) < 6){hp[i, t] <- hTobs[i, firstObs]} # Scenario 1

      nlast <- max(hh[i, 1:t] * tvec[1:t], na.rm = T)
      nnext <- min(hh[i, t:nt] * tvec[t:nt], na.rm = T)
      nlapsb <- t - nlast
      nlapsf <- nnext - t

      if(nlapsb <= 6 & nlapsf > 6){hp[i, t] <- hTobs[i, nlast]} # Scenario 3
      if(nlapsb > 6 & nlapsf <= 6){hp[i, t] <- hTobs[i, nnext]} # Scenario 4
      if(nlapsb > 6 & nlapsf > 6){next} # Scenario 5
#       if(t > lastObs & t - lastObs < 6){hp[i, t] <- hTobs[i, lastObs]} # Scenario 6

      if(nlapsb <= 6 & nlapsf <= 6){# Scenario 2
        tlast  <- hTobs[i, nlast ]
        tnext  <- hTobs[i, nnext ]
        h0     <- hTobs[i, nlast]
        if(!is.finite(nnext))nnext <- t
        if(!is.finite(tlast))tlast <- tnext
        if(!is.finite(tnext))tnext <- tlast
        if(!is.finite(nlast)){
          nlast <- t
          h0    <- tlast
        }

        slope  <- (tnext - tlast)/(nnext - nlast)
        hp[i, t] <- h0 + slope * (t - nlast)
      }
    }
  }

  round(hp,0)
}

initH <- function(hObs){

  hh    <- hObs*0 + 1
  tvec  <- 1:ncol(hObs)
  hTobs <- hObs
  hp    <- hTobs
  nt    <- ncol(hObs)
  n     <- nrow(hObs)

  for(i in 1:n){
    for(t in 1:nt){
      if(is.finite(hTobs[i, t]))next
      nlast <- max(hh[i, 1:t] * tvec[1:t], na.rm = T)
      nnext <- min(hh[i, t:nt] * tvec[t:nt], na.rm = T)
      tlast  <- hTobs[i, nlast ]
      tnext  <- hTobs[i, nnext ]
      h0     <- hTobs[i, nlast]
      if(!is.finite(nnext))nnext <- t
      if(!is.finite(tlast))tlast <- tnext
      if(!is.finite(tnext))tnext <- tlast
      if(!is.finite(nlast)){
        nlast <- t
        h0    <- tlast
      }

      slope  <- (tnext - tlast)/(nnext - nlast)
      hp[i, t] <- h0 + slope * (t - nlast)
    }
  }

  round(hp,0)
}

initH6months <- function(hObs){ # In contrast to initH(), this just interpolates if the sightings are within 6 months of a missing data value

  hh    <- hObs*0 + 1
  tvec  <- 1:ncol(hObs)
  hTobs <- hObs
  hp    <- hTobs
  nt    <- ncol(hObs)
  n     <- nrow(hObs)

  for(i in 1:n){
    firstObs <- min(which(is.finite(hTobs[i,])), na.rm = T)
    lastObs <- max(which(is.finite(hTobs[i,])), na.rm = T)

    for(t in 1:nt){
      if(is.finite(hTobs[i, t]))next

      if(t < firstSight[i]){next}
      if(t < firstObs & abs(firstObs - t) < 6){hp[i, t] <- hTobs[i, firstObs]} # Scenario 1

      nlast <- max(hh[i, 1:t] * tvec[1:t], na.rm = T)
      nnext <- min(hh[i, t:nt] * tvec[t:nt], na.rm = T)
      nlapsb <- t - nlast
      nlapsf <- nnext - t

      if(nlapsb <= 6 & nlapsf > 6){hp[i, t] <- hTobs[i, nlast]} # Scenario 3
      if(nlapsb > 6 & nlapsf <= 6){hp[i, t] <- hTobs[i, nnext]} # Scenario 4
      if(nlapsb > 6 & nlapsf > 6){next} # Scenario 5
      if(t > lastObs & t - lastObs < 6){hp[i, t] <- hTobs[i, lastObs]} # Scenario 6

      if(nlapsb <= 6 & nlapsf <= 6){# Scenario 2
      tlast  <- hTobs[i, nlast ]
      tnext  <- hTobs[i, nnext ]
      h0     <- hTobs[i, nlast]
      if(!is.finite(nnext))nnext <- t
      if(!is.finite(tlast))tlast <- tnext
      if(!is.finite(tnext))tnext <- tlast
      if(!is.finite(nlast)){
        nlast <- t
        h0    <- tlast
      }

      slope  <- (tnext - tlast)/(nnext - nlast)
      hp[i, t] <- h0 + slope * (t - nlast)
    }
   }
  }
  round(hp,0)
}

row2Mat <- function(vec){

  if(is.matrix(vec))return(vec)
  vn  <- names(vec)
  vec <- matrix(vec,1)
  colnames(vec) <- vn
  vec
}


myNewrmultinom <- function(size,p){

  #n multinomial r.v. for a n by ncol(p) matrix of probs
  #each row of p is a probability vector

  p <- row2Mat(p)

  n     <- nrow(p)
  J     <- ncol(p)

  if(length(size) == 1)size <- rep(size,n)

  jord  <- sample(J,J)    #randomize order

  p <- row2Mat(p[,jord])

  y <- yy  <- matrix(0,n,J)
  sizej <- size
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

  if(n == 1)y[,J] <- size - sum(y)
  if(n > 1) y[,J] <- size - rowSums(y)

  yy[,jord] <- y
  yy

}

#' Take the priors for regional movement and downweight them
#'
#' This function accepts two things: 1) a 9 * 9 matrix
#' of priors, and a 1 * 9 vector of (down) weights. The goal
#' of the function is to simply iterate over each column in the
#' matrix, which is a from column, and divide each element
#' in the column by the corresponding value in the vector
#'
#' @param xmat - the 9 * 9 array of priors for movement
#' @param yvec - a 1 * 9 vector of weights to divide
#'     the priors by
#' @return new - a 9 * 9 re-weighted movement prior
reWeight <- function(xmat, yvec){
   new <- t(apply(xmat, 1, function(x) x / yvec))
   new[new < 1] <- 1
   new
}
