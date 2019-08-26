if(RobHome) {
  nullfem <- read.csv(file = '/Users/robs/Documents/research/projects/PCAD/rightwhales/data/NulliparousForR.csv')
  delfem  <- read.csv(file = '/Users/robs/Documents/research/projects/PCAD/rightwhales/data/DelayedForR.csv')
  lgaps   <- read.csv(file = '/Users/robs/Documents/research/projects/PCAD/rightwhales/data/maxGapsforPhilip.csv')
  } else {
    nullfem <- read.csv(file = 'data/NulliparousForR.csv')
    delfem  <- read.csv(file = 'data/DelayedForR.csv')
    lgaps   <- read.csv(file = 'data/maxGapsforPhilip.csv')
  }



# load(file = 'allmerged.rdata') # Contains the merged distributions from the expert elicitation
if(subAn){load("/Users/rob/Documents/research/projects/PCAD/rightwhales/data/SubWhales.rdata")} # Contains the 50 whales to work on
load('egDemogData_2016-04-29.rdata')
# I'm going to remove 1005 from this as it's a MALE, but in calfTable
calves <- calves[!calves$EGNo == 1005,]
calfTable <- calves
deadTable <- read.csv('DeadYearMonthFeb2016.csv',header=T)
if(subAn){deadTable <- deadTable[deadTable$SightingEGNo %in% subWhales$EGNo,]}

sights <- data.frame(sights)
sights <- sights[is.finite(sights[, 'Date']), ]
# n.b. that here is the place we can take out some of the regions: X and maybe East
sights <- sights[!sights$RegionCode == 'X',]
sights <- sights[!sights$RegionCode == 'EAST',]
sights$RegionCode <- factor(sights$RegionCode)
#regNames <- regNames[which(regNames %in% unique(sights$RegionCode)) ]
if(subAn){sights <- sights[sights$SightingEGNo %in% subWhales$EGNo, ]}
if(subYrs){stopYr <- 2005} else {stopYr  <- max(sights$SightingYear, na.rm = T)}

if(remVHP){sights <- sights[as.Date(sights$Date) <= as.Date('2009-12-31'),] # remove the last year of sightings
          sights$BodyFatCode[sights$SightingYear == 2009] <- NA # remove an additional year of VHP data for each parameter
          sights$SkinCode[sights$SightingYear == 2009] <- NA
          sights$CyamidCode[sights$SightingYear == 2009] <- NA
          sights$LeftRakeMarkCode[sights$SightingYear == 2009] <- NA
          sights$RightRakeMarkCode[sights$SightingYear == 2009] <- NA
}

# figure out those alive before after 1970
bvec <- sights %>%
  group_by(SightingEGNo) %>%
  summarise(min_year = min(SightingYear, na.rm = T),
            max_year = max(SightingYear, na.rm = T)) %>%
  filter(min_year <= 1970 & max_year > 1970)


gender <- sights[, 'GenderCode']
# if(subYrs){stopYr <- 2005} else {stopYr  <- max(sights$SightingYear, na.rm = T)}
yrvec  <- startYr:stopYr
nyr    <- length(yrvec)

monYr  <- cbind( rep(c(1:12), nyr), rep(yrvec, each = 12) )
myName <- paste(monYr[, 1], monYr[, 2], sep = '-')
nt     <- nrow(monYr)
tIndex <- c(1:nt)

hseq <- c(0:100)

# load(file='pStates.rdata')
if(exists('breakEnt')){if(length(breakEnt == 2)){load('tStatesLumped.rdata')} else {load('tStates.rdata')}
                       m.idx <- which(colnames(tStates) %in% myName)
                       tStates <- tStates[, m.idx]}

# m.idx <- which(colnames(pStates) %in% myName)
# pStates <- pStates[, m.idx]

moveSeason <- rep(1,12)
moveSeason[3:4]  <- 2
moveSeason[5:6]  <- 3
moveSeason[7:10] <- 4

sightmy  <- paste(sights[, 'SightingMonth'], sights[, 'SightingYear'], sep = '-')
dataTime <- match(sightmy, myName)
sights   <- sights[is.finite(dataTime),] # Here's where they get pared down to the startYr
dataTime <- dataTime[is.finite(dataTime)]


dataNames <- c('BodyFatCode', 'LeftRakeMarkCode', 'RightRakeMarkCode',
               'SkinCode', 'CyamidCode')


dataMats <- numeric(0)
for(k in 1:length(dataNames)){

  dataMats <- append(dataMats, list(healthIndex(dataNames[k])))

}
dataNames <- c('BF', 'rakeL', 'rakeR', 'skin', 'cyam')
names(dataMats) <- dataNames

# next line is a test to see if otput is correct
# k <- 1
# CNAME = 'BodyFatCode'
# i <- 120
# wi <- which(sights[,'SightingEGNo'] == 1243)
# data.frame(orig = sights[wi, 'BodyFatCode'], interp = dataMats$BF[wi] )
# data.frame(orig = sights[wi, 'SkinCode'], interp = dataMats$skin[wi] )
# end test

# dataMats$rake <- dataMats$rakeL + dataMats$rakeR - 1
dataMats$rake <- apply(cbind(dataMats$rakeL, dataMats$rakeR), 1, function(x) mean(na.omit(x))) #- 1
dataMats$rake[which(dataMats$rake == 1.5)] <- 1
dataMats$rake[which(dataMats$rake == 2.5)] <- 2
dataNames <- c('BF', 'rake', 'skin', 'cyam')


stab   <- table(sights[, 'SightingEGNo'], sights[, 'RegionCode'],
                sights[, 'SightingYear'], sights[, 'SightingMonth'])


ID     <- sort(unique(sights[, 'SightingEGNo']))
n      <- length(ID)
regID  <- sort(unique(sights[, 'RegionCode']))
nreg   <- length(regID)

sightMat <- array(0, c(n, nt, nreg))
health <- matrix(NA, n, nt)
hStates <- health
rownames(health) <- ID
#colnames(health) <- paste(monYr[,1],monYr[,2],sep='-')
stage  <- health

rake <- skin <- cyam <- health

firstSight <- lastSight <- rep(NA, n)

calfByYr <- matrix(0, n, nyr)
colnames(calfByYr) <- yrvec
rownames(calfByYr) <- ID
calfByYr[cbind(match(calves[, 'EGNo'], ID), match(calves[, 'CalvingYear'], yrvec))] <- 1
wcy <- which(calfByYr == 1, arr.ind = T)

# creating two temporal indices (wcg, and wcr) to get the before calving year (wcg for gestation),
# and the after calving year (wcr, for recovery)
wcr <- wcg <- wcy
wcg <- wcg[wcg[, 2] != 1, ]
wcg[, 2] <- wcg[, 2] - 1
wcg[which(wcg[, 2] == 0), 2] <- 1

calfByYr[wcg] <- 2

calves  <- health * 0

#known deaths
deadTime <- match(paste(deadTable[, 'SightingMonth'], deadTable[, 'SightingYear'], sep = '-'), myName)
deadID   <- match(deadTable[, 'SightingEGNo'], ID)
deadReg  <- match(deadTable[, 'RegionCode'], regID)

gender <- rep(character(0), n)

for(i in 1:n){

  wi <- which(sights[, 'SightingEGNo'] == ID[i])

  gi <- sights[wi, 'GenderCode']

  gender[i] <- 'F'
  if('M' %in% gi)gender[i] <- 'M'

  stage[i, dataTime[wi]] <- sights[wi, 'AgeClassCode']

  hh  <- rep(NA, nt)
  hh[dataTime[wi]] <- dataMats$BF[wi]

  hStates[i, ] <- hh
  skin[i, dataTime[wi]]   <- dataMats$skin[wi]
  rake[i, dataTime[wi]]   <- dataMats$rake[wi]
  cyam[i, dataTime[wi]]   <- dataMats$cyam[wi]

  calfi <- calfByYr[i, ]
  w0 <- which(calfi == 1)
  yearC <- yrvec[w0]

  for(t in yearC){

    wt <- which(monYr[, 2] == t) # find the index of the months in the calving year
    wt <- c(wt[1]:(wt[1] + 9)) # this assigns 1s for the first 10 months of the calving year
    wtn <- c((wt[1] + 10):(wt[1] + 23))
    wtn <- wtn[wtn <= nt] # I put this in to deal with the sub animals that calve towards the end of the 1995-2005 period
    calves[i, wt] <- 1 # calving year
    calves[i, wtn] <- 3 # resting years

    if(t %in% startYr)next
    wt  <- which(monYr[, 2] == t) # find the index of the months in the calving year
    wtt <- c((wt[1] - 12):(wt[1] - 1)) # establish the pregnancy year
    calves[i, wtt] <- 2 # pregnancy year

  }

  # Put in the intervening months as 4's, i.e. available to calve
  if(length(yearC) > 0){
    c31 <- which(calves[i,] == 2)[1] # first month of the first pregnancy year
    c32 <- which(calves[i,] == 1)[length(which(calves[i, ] == 1))] # last month of the last calving year
    c3na <- which(is.na(calves[i, c31:c32])) # find any months between these that don't have a status
    calves[i, (c31:c32)[c3na]] <- 4 # put 4, or available, in these months
  }


  for(y in yrvec){
    ti <- which(monYr[,2] == y)
    yi <- which(dimnames(stab)[[3]] == y)
    if(length(yi) == 0)next

    for(r in 1:nreg){
      #     sightMat[i,ti,r] <- apply(stab[i,r,yi,,],1,sum)
      sightMat[i, ti, r] <- stab[i, r, yi, ]
    }
  }

  sm <- sightMat[i, , ]
  si <- apply(sm, 1, max)
  ws <- which(si > 0)
  t1 <- min(ws)
  t2 <- max(ws)
  w0 <- numeric(0)

  ss <- si
  ss[ss > 0] <- 1
  tmp <- c(1:length(si))*ss
  lastSight[i] <- max(tmp, na.rm = TRUE)
  tmp[tmp == 0] <- NA
  firstSight[i] <- min(tmp, na.rm = TRUE)
  if((lastSight[i] + missingInt) > nt){   wf <- firstSight[i]:nt} else {   wf <- firstSight[i]:(lastSight[i] + missingInt)}
  #    if(ID[i] %in% bvec$SightingEGNo){firstSight[i] <- 1}

#   # put the trailing 4's back in after the last calf, but before last sight
#   if(length(yearC) > 0){
#     c3max <- max(which(calves[i,] == 3))
#     c1max <- max(which(calves[i,] == 1))
#     calves[i, ((max(c3max, c1max)):lastSight[i])] <- 4
#   }

  if(entFlag){hEnt <- breaks2mids(breakEnt, tStates[i, wf])}
  if(pFlag){hProp  <- breaks2mids(breakProp, pStates[i, wf])}
  if(bfFlag){hBf   <- breaks2mids(breakBf, hh[wf])}
  if(skFlag){hSk   <- breaks2mids(breakSkinLims, skin[i, wf])}
  if(rkFlag){hRa   <- breaks2mids(breakRake, rake[i, wf])}
  if(cyFlag){hCy   <- breaks2mids(breakCyam, cyam[i, wf])}
  if(caFlag){hCa   <- breaks2mids(breakCalf, calves[i, wf])}

  mh   <- apply(cbind(if(exists('hEnt')){hEnt}, if(exists('hBf')){hBf}, if(exists('hSk')){hSk},
                      if(exists('hRa')){hRa}, if(exists('hCy')){hCy}, if(exists('hCa')){hCa}), 1, mean, na.rm = TRUE)
  mall <- rep(mean(mh, na.rm = TRUE), length(mh))

  wg <- which(is.finite(mh))
  if(length(wg) < 8){
    health[i, wf] <- mall
    next
  }

  ti <- c(1:nt)[wf]

  ytemp <- mh[is.finite(mh)]
  xtemp <- ti[is.finite(mh)]

  cc <- loess(ytemp ~  xtemp)
  mh <- predict(cc, data.frame(xtemp = ti))
  mall[is.finite(mh)] <- mh[is.finite(mh)]
  mh <- mall

  #  if(t2 < (nt - 20))mh[(t2+20):nt] <- NA
  #  mh[1:(t1-1)] <- NA

  health[i, wf] <- mh

}

# if remVHP, remove the last several years worth of VHP data to see its effect on the model
# if(remVHP)
# {
#   v.idx <- 301:432
#   hStates[, v.idx] <- NA
#   skin[, v.idx]    <- NA
#   rake[, v.idx]    <- NA
#   cyam[, v.idx]    <- NA
# }

# If subAn, Call a function to put in fake health data into the early parts of the trace
if(subAn & subAnImp)
{
  source('/Users/rob/Documents/code/rss10/rightwhales/subAnimalImpute.r')
}

#################################
# Missing Data Imputation Routine
#################################

notMissH <- which(is.finite(hStates))
missH <- which(is.na(hStates))
hobs <- hStates
hobs[missH] <- NA

notMissSkin <- which(is.finite(skin))
missSkin <- which(is.na(skin))
skobs <- skin
skobs[missSkin] <- NA

if(missDataImp){

  # Body Fat
  hTable <- rbind( c(5:1), c(1, 5, 10, 5, 1), c(1:5) )
  hTable[,3] <- c(1, 5, 1)
  hTable <- hTable / matrix( colSums(hTable), 3, 5 ,byrow = TRUE)
  ntab   <- ncol(hTable)

  hStates <- initH6months(hobs) # linearly interpolates on body fat scale
  hStatePrior <- makeHPrior6months(hobs, hStates, skinFlag = FALSE) # linearly interpolates on hTable scale
  pregidx <- which(calves == 2 & is.na(hobs) & hStatePrior != 5)
  hStatePrior[pregidx] <- 5

  # Skin
  skTable <- rbind( c(5:1), c(1:5) )
  skTable <- skTable / matrix( colSums(skTable), 2, 5 ,byrow = TRUE)

  ntab   <- ncol(skTable)

  skin <- initH(skobs) # linearly interpolates on skin scale
  skinStatePrior <- makeHPrior(skobs, skin, skinFlag = TRUE) # linearly interpolates on hTable scale

}

#####################################
# End Missing Data Imputation Routine
#####################################

sightAll <- apply(sightMat, 1:2, sum)  #individual sightings by time



iByr <- matrix(c(1:nreg), n, nreg, byrow = TRUE)  #move this


effort <- data.frame(effort)

effortmy  <- paste(effort[, 'month'], effort[, 'year'], sep = '-')
emat      <- effort[effortmy %in% myName, ]
effortmy  <- effortmy[effortmy %in% myName]

regNum   <- match(emat[, 'region'], regID)
regNum[is.na(regNum)] <- 1                    #temporary!

effortmat <- matrix(0, nt, nreg)
colnames(effortmat) <- regID
rownames(effortmat) <- myName
eIndex    <- cbind(match(effortmy, myName), regNum)

effortmat[eIndex] <- emat[, 'effort']


effortState <- sightMat * 0
for(i in 1:n)effortState[i, , ] <- effortmat

# mean(sightMat[effortState > 0]/effortState[effortState > 0],na.rm=T)


# Get the max gap from lgaps ordered along with ID & catch the cases (NA and -Inf)
# For those cases, set to the median of all animals
iddat <- merge(as.data.frame(ID), lgaps, by.x = 'ID', by.y = 'SightingEGNo', all.x = TRUE)
iddat[which(!is.finite(iddat$maxgap)), 'maxgap'] <- median(iddat$maxgap, na.rm = T)
iddat[which(is.na(iddat$maxgap)), 'maxgap'] <- median(iddat$maxgap, na.rm = T)

