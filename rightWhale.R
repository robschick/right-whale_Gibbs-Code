# this is the version in the git archive
rm(list=ls())
RobHome <- FALSE
JimHome <- F
library(compiler)
library(matrixStats)
library(stringr)
library(eliciteg)
library(dplyr)
options(error = dump.frames)
subAn        <- FALSE # If TRUE, just 50 animals
subAnImp     <- FALSE # If TRUE, add in fake data to inspect early parts of health traces
subYrs       <- FALSE # If TRUE, use just 1995-2005
FixSurv      <- FALSE # If TRUE, fix survival it at 1
Survhack     <- FALSE # If TRUE, use the NEAq presumed dead effect
REGINTERCEPT <- FALSE # If TRUE, include an intercept for all regions in the survival estimation
missDataImp  <- TRUE  # if TRUE the impute missing data for Body Fat and Skin
movePsa      <- FALSE # if TRUE, run a movement prior sensitivity analysis
remVHP       <- FALSE # if TRUE, remove the VHPs from 1995 to 2005
flatPrior    <- FALSE # if TRUE, every value in the Dirichlet is 1
longGap      <- FALSE # if TRUE, missingInt is == longest sighting gap of each individual
downWeight   <- FALSE  # if TRUE, use Philip Hamilton's proposed downweights to create new movement priors
multExp      <- FALSE #if true, use the multiple expert sampling
iter <- 1 # to allow me to only save in the second iteration, i.e. when BIG is running.
if(movePsa){movePsaVal <- 2}
if(RobHome) {

  wkdir <- '/Users/robs/Documents/research/projects/PCAD/rightwhales'
  setwd(wkdir)
  if(subYrs){startYr <- 1995} else {startYr <- 1970}
  source('/Users/robs/Documents/research/code/rightwhales/clarkFunctions.R')
  source('/Users/robs/Documents/research/code/rightwhales/rightWhaleFunctions.R')

} else {

  wkdir <- '/home/rob/rightwhales'
  setwd(wkdir)
  if(subYrs){startYr <- 1995} else {startYr <- 1970}
  source('clarkFunctions.R')
  source('rightWhaleFunctions.R')

}

modelname <- 'eg_213_'
mname <- new.env()
assign('modelname', value = modelname, envir = mname)
ng <- 50000
mname$modelname <- paste(mname$modelname, 'ng_', ng, sep='')
if(file.exists(paste("gibbsoutput/", mname$modelname, '.rdata', sep="") ) ){stop('Output from a model with that name already exists')}
bfFlag  <- TRUE
entFlag <- FALSE
rkFlag  <- TRUE
skFlag  <- TRUE
cyFlag  <- TRUE
caFlag  <- FALSE
pFlag   <- FALSE
hnamevec <- c('fat', 'ent', 'rake', 'skin', 'cyam', 'calf', 'prop')
truevec  <- c(bfFlag, entFlag, rkFlag, skFlag, cyFlag, caFlag, pFlag)
hnamevec <- hnamevec[truevec]

missingInt <- 48 # time (in months) for which the animal can still be alive after last sighting


hoffset    <- 5 # factor by which health can change in a month
JUV        <- TRUE
if(JUV){juvAge <- 9} else {juvAge <- 0}  # age when juveniles become adults
bnt        <- 100 # time offset to give to animals that are still alive

#body fat
if(bfFlag){
  lo <- c(20, 70)
  hi <- c(40, 100)
  breakBfLims <- rbind(lo, hi)
  breakBf     <- apply(breakBfLims, 2, mean)
}

#entanglement
if(entFlag){
  lo <- c(10, 50)
  hi <- c(80, 90)
  breakEntLims <- rbind(lo, hi)
  breakEnt <- apply(breakEntLims, 2, mean)
}

#propeller
if(pFlag){
  lo <- c(10, 30, 60)
  hi <- c(80, 90, 100)
  breakPropLims <- rbind(lo, hi)
  breakProp <- apply(breakPropLims, 2, mean)
}

#rakes
if(rkFlag){
  lo <- c(23, 55)
  hi <- c(75, 90)
  breakRakeLims <- rbind(lo, hi)
  breakRake <- apply(breakRakeLims, 2, mean)
}

#skin
if(skFlag){
  lo <- 30
  hi <- 100
  breakSkinLims <- matrix(c(lo, hi), 2, 1); rownames(breakSkinLims) <- c('lo', 'hi')
  breakSkin <- mean(breakSkinLims)
}

#cyam
if(cyFlag){
  lo <- 20
  hi <- 80
  breakCyamLims <- matrix(c(lo, hi), 2, 1); rownames(breakCyamLims) <- c('lo', 'hi')
  breakCyam <- mean(breakCyamLims)
}

#calf
if(caFlag){
  lo <- 25
  hi <- 85
  breakCalfLims <- rbind(lo, hi)
  breakCalf <- apply(breakCalfLims, 2, mean)
}

if(RobHome){source('/Users/robs/Documents/research/code/rightwhales/rightWhaleInput.R')} else {source('rightWhaleInput.R')}



#survival variables
varS <- c('health',regID)

pint <- -2
h50 <- 40   #value of health p(surv) = .5
maxSurv <- .99 #maximum survival rate when health = 100
ph <- -pint/h50
ps <- length(varS)
priorBSurv  <- matrix(c(ph, rep(pint, nreg)), ncol = 1)
rownames(priorBSurv) <- varS
priorVS <- diag(1, ps)
#loBS    <- c(.3,rep(-4.5,nreg))
#hiBS    <- c(100,rep(6,nreg))

bg      <- priorBSurv

pv1 <- log(.002)
pv2 <- log(.0005)


#health state variable parameters
varQ <- c('intercept', 'h_t-1', 'juv', 'age', 'age^2')
pq   <- length(varQ)
priorBQ <- matrix(c(0, 1, -.1, .1, -.0005), pq, 1)
rownames(priorBQ) <- varQ
priorIVQ <- solve(diag(c(10, 1/n/nt, 10, 10, 10)))       #strong prior is low value
loBQ    <- c(-6, -2, -.5, 0, -6)                          #for truncated prior
hiBQ    <- c(6, 2, 1.5, 6, 0)
bq      <- priorBQ

if(!JUV){
  priorBQ[varQ == 'juv'] <- 0
  priorIVQ[c(varQ == 'juv', varQ == 'juv')] <- 0
  bq <- priorBQ
  loBQ[varQ == 'juv'] <- 0
  hiBQ[varQ == 'juv'] <- 0
}


muh    <- 10       #mean process variance for health, monthly scale
h1     <- nt*n
h2     <- muh*(h1 - 1)
herror <- muh

#detection: mean sighting = n/no. regions/effort
lamMean <- sum(sightMat,na.rm=T)/n/sum(effortmat,na.rm=T) * 100  #weighting
aLam <- 3
bLam <- aLam*(1/lamMean - 1)

bLam1  <- bLam
amean <- bLam*lamMean
aLam1 <- bLam1*amean

lamda <- rep(lamMean,n)

# ###########################
# Updated (12/16/11) movements

# load in new adjusted priors (if needed). (Note that implementation of the
# new priors means two things: 1) using the new adjusted values from Philip,
# and 2) using the down weights, which will come in below.)

if(downWeight) {
  adultMalesNew <- read.csv(file = 'data/adultMales_PhilipUpdated_Sep2016.csv', header = T)
  adultFemalesNew <- read.csv(file = 'data/adultFemales_PhilipUpdated_Sep2016.csv', header = T)
  xGenderNew <- read.csv(file = 'data/xGender_PhilipUpdated_Sep2016.csv', header = T)
  juvMalesNew <- read.csv(file = 'data/juvMales_PhilipUpdated_Sep2016.csv', header = T)
  juvFemalesNew <- read.csv(file = 'data/juvFemales_PhilipUpdated_Sep2016.csv', header = T)
}

wt <- n    # weight for prior; n is the baseline
if(movePsa){wt <- n * movePsaVal}
if(downWeight) {
  moveAll <- InitialMove(adultMalesNew, wt)
  } else {
  moveAll <- InitialMove(adultMales, wt)
}
maleMoveProb <- moveAll$moveProb
maleMove     <- moveAll$move
if(flatPrior){maleMove[maleMove > 0] <- 1}

if(downWeight) {
  moveAll <- InitialMove(adultFemalesNew, wt) # using new adjust priors
} else {
  moveAll <- InitialMove(adultFemales, wt) # using old priors
}
femMoveProb <- moveAll$moveProb
femMove     <- moveAll$move
if(flatPrior){femMove[femMove > 0] <- 1}

if(downWeight) {
  moveAll <- InitialMove(xGenderNew, wt)
} else {
  moveAll <- InitialMove(xGender, wt)
}
xMoveProb <- moveAll$moveProb
xMove     <- moveAll$move
if(flatPrior){xMove[xMove > 0] <- 1}

if(downWeight) {
  moveAll <- InitialMove(juvFemalesNew, wt)
} else {
  moveAll <- InitialMove(juvFemales, wt)
}
juvMoveProb <- moveAll$moveProb
juvMove     <- moveAll$move
if(flatPrior){juvMove[juvMove > 0] <- 1}

# Downweighting the priors occurs here. That is we use Philip Hamilton's (NEAq)
# revised priors, both in terms on individual values within each monthly array,
# and the month/gender/region specific down-weighting that he proposed:

if(downWeight){

  if(RobHome) {
    pw <- read_csv('/Users/robs/Documents/research/projects/midaMove/data/PhilipsMay2015revisedPriorWeightsLongFormat.csv')
  } else {
    pw <- read_csv('data/PhilipsMay2015revisedPriorWeightsLongFormat.csv')
  }

  for(i in seq_along(1:dim(maleMove)[3])){
    df <- maleMove[, , i]
    wvec <- filter(pw, gender == 'adMale' & monthNum == i) %>% .$DownWeight
    maleMove[, , i] <- reWeight(df, wvec)
  }

  for(i in seq_along(1:dim(femMove)[3])){
    df <- femMove[, , i]
    wvec <- filter(pw, gender == 'adFemale' & monthNum == i) %>% .$DownWeight
    femMove[, , i] <- reWeight(df, wvec)
  }

  for(i in seq_along(1:dim(juvMove)[3])){
    df <- juvMove[, , i]
    wvec <- filter(pw, gender == 'juveniles' & monthNum == i) %>% .$DownWeight
    juvMove[, , i] <- reWeight(df, wvec)
  }

}

par(mfrow=c(3,2))
#health imputation


varInt <- 10    #variance for intercepts
varSlp <- 1/n   #variance for slopes

# body fat
if(bfFlag){
  tmp         <- makePriorB(c(10, 30), c(-1, -.0001), breakBf, varInt = varInt, varSlp = varSlp)
  priorBH     <- tmp$meanPrior
  priorVBH    <- tmp$varPrior
  bcovHealth  <- diag(.01,length(breakBf)+1)
  bh    <- priorBH
  plotLogit(breakBfLims,breakBf,hseq,priorBH)
  title('Body Condition')
  #   dev.off()
}

#entanglement
if(entFlag){
  tmp    <- makePriorB(c(10,30), c(-.5,-.2), breakEnt, varInt = varInt, varSlp = varSlp)
  priorBE <- tmp$meanPrior
  priorVBE <- tmp$varPrior
  bcovEnt  <- diag(.0001,length(breakEnt)+1)
  bEnt  <- priorBE
  plotLogit(breakEntLims,breakEnt,hseq,priorBE)
  title('entanglement')
}

# Propeller
if(pFlag){
  tmp    <- makePriorB(c(10,30), c(-.5,-.01), breakProp, varInt = varInt, varSlp = varSlp)
  priorBP <- tmp$meanPrior
  priorVBP <- tmp$varPrior
  bcovProp  <- diag(.0001,length(breakProp)+1)
  bProp  <- priorBP
  plotLogit(breakPropLims, breakProp, hseq, priorBP)
  title('Prop Scarring')
}

#rake scale
if(rkFlag){
  tmp  <- makePriorB(c(10, 30), c(-1, -.4), breakRake, varInt = varInt, varSlp = varSlp)
  priorBR <- tmp$meanPrior
  priorVBR <- tmp$varPrior
  bcovRake <- diag(.0000001,(length(breakRake)+1))
  bRake <- priorBR
  plotLogit(breakRakeLims,breakRake,hseq,priorBR)
  title('rake injury')
}

#skin scale
if(skFlag){
  bcovSkin <- diag(c(.01,.001))
  b1 <- -.25
  b0 <- -b1*breakSkin
  priorBS <- matrix(c(b0,b1),1,2)
  priorVBS <- diag(1/n,2)
  bSkin <- priorBS
  plotLogit(breakSkinLims,breakSkin,hseq,priorBS)
  title('skin injury')
}

#cyam scale
if(cyFlag){
  bcovCyam <- diag(c(.001,.0001))
  b1 <- -.25
  b0 <- -b1*breakCyam
  priorBC <- matrix(c(b0,b1),1,2)
  priorVBC <- diag(1/n,2)
  bCyam <- priorBC
  plotLogit(breakCyamLims,breakCyam,hseq,priorBC)
  title('cyamids')
}

#calves
if(caFlag){
  tmp   <- makePriorB(c(10, 30), c(-.25, -.15), breakCalf, varInt = varInt, varSlp = varSlp)
  priorBF <- tmp$meanPrior
  priorVBF <- tmp$varPrior
  bcovCalf <- diag(.01,1)  #break[1:2] and steepness
  bCalf <- priorBF
  plotLogit(breakCalfLims,breakCalf,hseq,priorBF)
  title('calving status')
}

tiny <- 1e-20

#initialize latent states
rsState <- sightMat
rsState[sightMat > 0] <- 1  #known alive and location
rState <- rsState           #currently imputed location
firstYr <- rep(NA,n)

xinit <- yinit <- numeric(0)
death <- rep(nt + bnt, n)


for(i in 1:n){

  si <- apply(sightMat[i,,],1,max)
  ws <- which(si > 0)
  t1 <- min(ws)
  t2 <- max(ws)

  firstYr[i] <- monYr[t1,2]
  wdead <- which(deadID == i)

  rState[i,t1,rsState[i,t1,] == 1] <- 1

  for(t in t1:nt){

    if(t > nt)break
    surv0 <- 1
    surv1 <- 1
    age <- (t - t1)/12
    ja  <- 1
    if(age > juvAge)ja <- 0

    health1 <- health[i, t-1]
    if(subAn){if(t == 1){health1 <- 30 * mean(hStates[i,] , na.rm=T)}}
    if(is.na(health1)) health1 <- 30 * mean(hStates[i,] , na.rm=T)
    if(is.na(health1)) health1 <- 75 # to error trap for animals with no body condition scores, e.g. 1002


    health[i,t] <- c(1, health1, ja, age, age ^ 2) %*% bq
    if(health[i, t] > 90)health[i, t] <- 90
    if(health[i, t] < 10)health[i, t] <- 10

    X <- c(health1, rState[i, t, ])

    xinit <- rbind(xinit, X)
    yinit <- c(yinit, 1)

    if(t > t2){

    surv0 <- inv.logit(X%*%priorBSurv)
    if(is.na(surv0))surv0 <- 1
    survive <- rbinom(1,1,surv0)
    if(length(wdead) > 0){
    if(t <  deadTime[wdead])survive <- 1
    if(t >= deadTime[wdead])survive <- 0
    }

    if(survive == 0 | t > (lastSight[i] + missingInt)){
    death[i] <- t
    yinit[length(yinit)] <- 0
    health[i,t] <- 0
    break
    }
    }

    w0 <- which(sightMat[i,t,] > 0)
    if(length(w0) > 0){
      rsState[i,t,w0] <- 1
      xinit[nrow(xinit),(w0+1)] <- 1
      next
    }

    w1 <- which(rState[i,t-1,] == 1)

    monNow <- moveSeason[monYr[t,1]]

    pt <- rep(1,nreg)
    if(length(w1) > 0){
      pt <- femMove[,w1[1],monNow]
      if(gender[i] == 'M')pt <- maleMove[,w1[1],monNow]
      pt <- pt*exp(-lamda[i]*effortmat[t,])
    }
    if(t <= t2){
      w2 <- which(rState[i,t+1,] == 1)
      if(length(w2) > 0){
        ps <- femMove[w2[1],,monNow]*pt
        if(gender[i] == 'M')ps <- maleMove[w2[1],,monNow]*pt
        pt <- ps
      }
      #       if(length(w2) > 0)pt <- moveProb[w2[1],]*pt
    }
    pt <- pt/sum(pt)
    q  <- which(rmultinom(1,1,pt) == 1)
    rState[i,t,q] <- 1

    xinit[nrow(xinit),(q+1)] <- 1
    # print(c(t, death[i]))
  }
  #    print(paste(ID[i], round(health[i,t1], 1), round(mean(hStates[i,] , na.rm=T), 2), sep = ' : '))
  #  print(c(t, death[i]))
}

#how often is an individual sighted in more than one region within a month?

tmp <- sightMat
tmp[tmp > 0] <- 1
multSite <- rep(NA,n)
for(i in 1:n){
  multSite[i] <- length(which(apply(tmp[i,,],1,sum) > 1))
}


pvals <- c(0.1, 0.99)              #2 probabilities for 2 health values
hvals <- c(15, 90)
pvals <- log(pvals/(1 - pvals))
b0    <- (pvals[1] * hvals[2] - pvals[2] * hvals[1]) / diff(hvals)
b1    <- diff(pvals) / diff(hvals)

bg <- matrix(c(b1, rep(b0, nreg)), ncol = 1)
loBS <- bg * 1.2
hiBS <- bg * 0.3

loBS[1] <- 0.8 * b1
hiBS[1] <- 0.2
hiBS[2:nrow(bg)] <- 0

propBS <- .001*solve(crossprod(xinit))

#####################################
ng <- ng
if(RobHome){source('/Users/robs/Documents/research/code/rightwhales/writeModelSummary.R')} else {source('writeModelSummary.R')}


if(!subAn){death[deadID] <- deadTime}

if (longGap) {
  missingInt <- iddat$maxgap
  death[death > (lastSight + missingInt)] <- lastSight[death > (lastSight + missingInt)] + missingInt[death > (lastSight + missingInt)]
} else {
  death[death > (lastSight + missingInt)] <- lastSight[death > (lastSight + missingInt)] + missingInt
}


save.image(file = paste('UpThroughInitLoop_', modelname,
                        noquote(as.character(Sys.Date())), '.rdata', sep = '') )
if(RobHome) {

  if(subAn){
    source('/Users/robs/Documents/research/code/rightwhales/rightWhaleGibbsSubAn.r')
  } else {
    source('/Users/robs/Documents/research/code/rightwhales/rightWhaleGibbs.R')
    }
}

if(!RobHome) {
  source('rightWhaleGibbs.r')
}


ngg <- g-1
save.image(file=paste("gibbsoutput/", mname$modelname, '.rdata', sep=""))


# # Run the next bit of the loop post-burnin
ng <- 50000
iter <- 2
mname$modelname <- paste(mname$modelname, '_BIG_', ng, sep = '')
if(RobHome) {

  if(subAn){
    source('/Users/robs/Documents/research/code/rightwhales/rightWhaleGibbsSubAn.r')
  } else {
    source('/Users/robs/Documents/research/code/rightwhales/rightWhaleGibbsLoop2.R')
    }
}
if(!RobHome) {
  source('rightWhaleGibbs.r')
}
ngg <- g-1
save.image(file = paste("gibbsoutput/", mname$modelname, '.rdata', sep = ""))
