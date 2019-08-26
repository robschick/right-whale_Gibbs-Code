#  # source('rightWhaleFunctions.r')
# # time1 <- Sys.time()
# hseq <- c(0:100)
# hseq.idx <- seq(1000, (ng - 30), by = 38) #every 38th one gives me 5k output

rsum <- rsState*0
deathyr <- matrix(0, n, (nt + bnt))
rownames(deathyr) <- ID
suml <- rep(0, n)
summ <- sumf <- sumx <- summ2 <- sumf2 <- sumx2 <- maleMoveProb*0
sumh <- matrix(0, n, nt)
sumh2 <- sumh
# Old Population Delineation by sub-category
# newhgibbsFemales <- newhgibbsMales <- newhgibbsJuvs <- newhgibbsRepFem <- newhgibbsRepFemBOF <- newhgibbsRepFemNonBOF <- sumh
# newhgibbsJuvs2 <- newhgibbsRepFem2 <- newhgibbsRepFemBOF2 <- newhgibbsRepFemNonBOF2 <- sumh2
#
# # Definition of Above Placeholders
# # newhgibbsMales # Males
# # newhgibbsJuvs # Juveniles (2-9)
# # newhgibbsRepFem # Reproductively Active Females
# # newhgibbsRepFemBOF # 'regular' Females, i.e. Go to BOF
# # newhgibbsRepFemNonBOF # irregular Females, i.e. Do not go to BOF

# New version of above delineation
#  Juveniles (1-2 years old)
#  Juveniles (3-9 years old)
#  Adult Males (>9 years old)
#  Adult Nulliparous Females
#  Adult Females with Delayed First Calf
#  Pregnant (year before calving)
#  Lactating (calving year)
#  Available to calve (Resting plus one or more years)
#  Resting November of lactating year through December of following year (14 months)

newhgibbsJuvs1 <- newhgibbsJuvs3 <- newhgibbsMales <- newhgibbsAdMales <- newhgibbsNonRepFem <- newhgibbsDelRepFem <- newhgibbsPregFem <- newhgibbsLacFem <- newhgibbsAvCavFem <- newhgibbsRestFem <- sumh
newhgibbsJuvs12 <- newhgibbsJuvs32 <- newhgibbsMales2 <- newhgibbsAdMales2 <- newhgibbsNonRepFem2<- newhgibbsDelRepFem2 <- newhgibbsPregFem2 <- newhgibbsLacFem2 <- newhgibbsAvCavFem2 <- newhgibbsRestFem2 <- sumh2

sumLamda <- sumLamda2 <- rep(0,n)

accepts <-  0

bggibbs   <- matrix(NA,ng,length(varS)); colnames(bggibbs) <- varS
bqgibbs   <- matrix(NA,ng,length(varQ)); colnames(bqgibbs) <- varQ
hhgibbs   <- matrix(NA,ng,2); colnames(hhgibbs) <- c('break1','break2')

lgibbs <- rep(NA,ng)
agibbs <- rep(NA,ng)
if(bfFlag){bhgibbs  <- matrix(NA,ng,4); colnames(bhgibbs) <- c('c01','c02','c11','c12')}
if(rkFlag){brgibbs  <- matrix(NA,ng,length(bRake))}
if(skFlag){bsgibbs  <- matrix(NA,ng,length(bSkin))}
if(cyFlag){bcgibbs  <- matrix(NA,ng,length(bCyam))}
if(caFlag){bfgibbs  <- matrix(NA,ng,length(bCalf))}
if(entFlag){begibbs <- matrix(NA,ng,length(bEnt))}

pFlag <- FALSE
if(pFlag){bpgibbs <- matrix(NA,ng,length(bProp))}

pv1 <- pv2 <- log(.01)

accept <- rep(0, 7)
names(accept) <- c('bf','rake','skin','cyam','calf','ent', 'prop')
htime <- matrix(1:nt,n,nt,byrow=T)


# Progress Bar
pb <- txtProgressBar(min = 0, max = ng, style = 3)

## Getting the Sub-groups health tabulated
# # 1) Males
# m.idx <- ddply(sights[sights$GenderCode == 'M',], .(SightingEGNo),
#                summarise, min_date = min(Date, na.rm = T),
#                max_date = max(Date, na.rm = T)) # Get the animals and times
#
# m.idx$MinDate <- paste(ifelse(format(m.idx$min_date, '%m') < 10,
#                               str_sub(format(m.idx$min_date, '%m'),2),
#                               format(m.idx$min_date, '%m')),
#                        format(m.idx$min_date, '%Y'), sep = '-')
#
# m.idx$MaxDate <- paste(ifelse(format(m.idx$max_date, '%m') < 10,
#                               str_sub(format(m.idx$max_date, '%m'),2),
#                               format(m.idx$max_date, '%m')),
#                        format(m.idx$max_date, '%Y'), sep = '-')
#
# m.idx$MinDateInt <- match(m.idx$MinDate, myName)
# m.idx$MaxDateInt <- match(m.idx$MaxDate, myName)
#
#
# # 2) Juvenile Females
# jd1.idx <- ddply(sights[sights$AgeClassCode == 'J' & sights$GenderCode == 'F',], .(SightingEGNo),
#                  summarise, min_date = min(Date, na.rm = T),
#                  max_date = max(Date, na.rm = T)) # Get the animals and times
# # Get the times in the right format
# jd1.idx$MinDate <- paste(ifelse(format(jd1.idx$min_date, '%m') < 10,
#                                 str_sub(format(jd1.idx$min_date, '%m'),2),
#                                 format(jd1.idx$min_date, '%m')),
#                          format(jd1.idx$min_date, '%Y'), sep = '-')
# jd1.idx$MaxDate <- paste(ifelse(format(jd1.idx$max_date, '%m') < 10,
#                                 str_sub(format(jd1.idx$max_date, '%m'),2),
#                                 format(jd1.idx$max_date, '%m')),
#                          format(jd1.idx$max_date, '%Y'), sep = '-')
#
# jd1.idx$MinDateInt <- match(jd1.idx$MinDate, myName)
# jd1.idx$MaxDateInt <- match(jd1.idx$MaxDate, myName)
#
# m.idx <- rbind(m.idx, jd1.idx) # bind them together to represent the 'Population level health'

# 3) UPDATE - 17 March 2014 - This is an almagamation of what I had previously. This is to get three groups anew: Adult Males, Young Juveniles, & Old Juveniles. I'm cutting and pasting this from a github commit SHA 826d3ed3 on the 3rd of February 2104
# 3.1) Adult Males
# am.idx <- ddply(sights[sights$AgeClassCode == 'A' & sights$GenderCode == 'M',], .(SightingEGNo),
#                summarise, min_date = min(Date, na.rm = T),
#                max_date = max(Date, na.rm = T)) # Get the animals and times
# Update 14 December 2015 - changing from plyr syntax to dplyr
am.idx <- sights %>%
  filter(AgeClassCode == 'A' & GenderCode == 'M') %>%
  group_by(SightingEGNo) %>%
  summarise(min_date = min(Date, na.rm = T),
            max_date = max(Date, na.rm = T))

am.idx$MinDate <- paste(ifelse(format(am.idx$min_date, '%m') < 10,
                              str_sub(format(am.idx$min_date, '%m'),2),
                              format(am.idx$min_date, '%m')),
                       format(am.idx$min_date, '%Y'), sep = '-')

am.idx$MaxDate <- paste(ifelse(format(am.idx$max_date, '%m') < 10,
                              str_sub(format(am.idx$max_date, '%m'),2),
                              format(am.idx$max_date, '%m')),
                       format(am.idx$max_date, '%Y'), sep = '-')

am.idx$MinDateInt <- match(am.idx$MinDate, myName)
am.idx$MaxDateInt <- match(am.idx$MaxDate, myName)

# 3.2) Young Juveniles
# jd1.idx <- ddply(sights[sights$AgeClassCode == 'J' & as.numeric(levels(sights$Age))[sights$Age] <= 2,], .(SightingEGNo),
#                  summarise, min_date = min(Date, na.rm = T),
#                  max_date = max(Date, na.rm = T)) # Get the animals and times

jd1.idx <- sights %>%
  filter(AgeClassCode == 'J' & as.numeric(Age) <= 2) %>%
  group_by(SightingEGNo) %>%
  summarise(min_date = min(Date, na.rm = T),
            max_date = max(Date, na.rm = T))

# Get the times in the right format
jd1.idx$MinDate <- paste(ifelse(format(jd1.idx$min_date, '%m') < 10,
                                str_sub(format(jd1.idx$min_date, '%m'),2),
                                format(jd1.idx$min_date, '%m')),
                         format(jd1.idx$min_date, '%Y'), sep = '-')
jd1.idx$MaxDate <- paste(ifelse(format(jd1.idx$max_date, '%m') < 10,
                                str_sub(format(jd1.idx$max_date, '%m'),2),
                                format(jd1.idx$max_date, '%m')),
                         format(jd1.idx$max_date, '%Y'), sep = '-')

jd1.idx$MinDateInt <- match(jd1.idx$MinDate, myName)
jd1.idx$MaxDateInt <- match(jd1.idx$MaxDate, myName)

# 3.3) Old Juveniles
# jd2.idx <- ddply(sights[sights$AgeClassCode == 'J'  & as.numeric(levels(sights$Age))[sights$Age] > 2,], .(SightingEGNo),
#                  summarise, min_date = min(Date, na.rm = T),
#                  max_date = max(Date, na.rm = T)) # Get the animals and times
jd2.idx <- sights %>%
  filter(AgeClassCode == 'J' & as.numeric(Age) > 2) %>%
  group_by(SightingEGNo) %>%
  summarise(min_date = min(Date, na.rm = T),
            max_date = max(Date, na.rm = T))

# Get the times in the right format
jd2.idx$MinDate <- paste(ifelse(format(jd2.idx$min_date, '%m') < 10,
                                str_sub(format(jd2.idx$min_date, '%m'),2),
                                format(jd2.idx$min_date, '%m')),
                         format(jd2.idx$min_date, '%Y'), sep = '-')
jd2.idx$MaxDate <- paste(ifelse(format(jd2.idx$max_date, '%m') < 10,
                                str_sub(format(jd2.idx$max_date, '%m'),2),
                                format(jd2.idx$max_date, '%m')),
                         format(jd2.idx$max_date, '%Y'), sep = '-')

jd2.idx$MinDateInt <- match(jd2.idx$MinDate, myName)
jd2.idx$MaxDateInt <- match(jd2.idx$MaxDate, myName)

m.idx <- rbind(am.idx, jd2.idx) # bind Adult males and Old Juveniles together to represent the 'Population level health'

# 4)  Adult Nulliparous Females
nullfem$MinDate    <- paste(1, nullfem$MinYear, sep = '-')
nullfem$MaxDate    <- paste(12, nullfem$MaxYear, sep = '-')
nullfem$MinDateInt <- match(nullfem$MinDate, myName)
nullfem$MaxDateInt <- match(nullfem$MaxDate, myName)
nullfem$MaxDateInt[is.na(nullfem$MaxDateInt)] <- nt

# 5) Adult Females with Delayed First Calf
# delfem            <- delfem[which(delfem$DiffFlag > 0),]
delfem$MinDate    <- paste(1, delfem$MinYear, sep = '-')
delfem$MaxDate    <- paste(12, delfem$MaxYear, sep = '-')
delfem$MinDateInt <- match(delfem$MinDate, myName)
delfem$MaxDateInt <- match(delfem$MaxDate, myName)
delfem$MaxDateInt[is.na(delfem$MaxDateInt)] <- nt

# 6) Pregnant Females - the logic here is to find all the times NOT pregnant, or for #7, 8, etc., lactating, available, etc. and set those to NA below in the Gibbs loop
fp.idx <- which(calves != 2 | is.na(calves), arr.ind = TRUE)

# 7) Lactating females
fl.idx <- which(calves != 1 | is.na(calves), arr.ind = TRUE)

# 8) Available to Calve
fac.idx <- which(calves != 4 | is.na(calves), arr.ind = TRUE)

# 9) Resting
fr.idx <- which(calves != 3 | is.na(calves), arr.ind = TRUE)

dev <- 0

for(g in 1:ng){
  setTxtProgressBar(pb, g)

  #update health before health parameters
  tmp <- updateHealth()
  health  <- tmp$health
  hStates <- tmp$hStates
  skin    <- tmp$skin
  acceptH <- tmp$accept

  tmp <- healthStateSpacePars()
  bq <- tmp$bq
  herror <- tmp$herror         #process error for health

  tmp <- updateStates()
  rState <- tmp$rs
  death   <- tmp$death
  deathyr[cbind(c(1:n),death)] <- deathyr[cbind(c(1:n),death)] + 1
  rsum <- rsum + rState
  maleMoveProb <- tmp$mm
  femMoveProb  <- tmp$mf
  xMoveProb    <- tmp$mx
  juvMoveProb  <- tmp$mj
  bg           <- tmp$bg
  #   accepts      <- accepts + tmp$accept

  if(g %in% c(50,100,200,400,600)){
    sss <- cov(bggibbs[1:(g-1),])
    testv <- try(solve(sss),T)
    if(!inherits(testv,'try-error'))propBS <- sss*.1
  }

  likeg <- 0

  # for(k in 1:3){
  #   par(mfrow=c(3,2))
  if(bfFlag){
    tmp <- healthPars(bh, priorBH, priorVBH, breakBf, breakBfLims, hobs, health, accept['bf'])
    bh  <- tmp$bgg
    breakBf     <- tmp$breakg
    likeg <- likeg + tmp$loglik
    # plotLogit(breakBfLims,breakBf,hseq,bh)
    aa     <- tmp$accept
    if(aa == 1)accept['bf'] <- 0
    if(aa == 0)accept['bf'] <- accept['bf'] + 1
  }

#   if(entFlag){
#     tmp <- healthPars(bEnt, priorBE, priorVBE, breakEnt, breakEntLims, tStates, health, accept['ent'])
#     bEnt      <- tmp$bgg
#     breakEnt <- tmp$breakg
#     likeg <- likeg + tmp$loglik
#     #  plotLogit(breakEntLims,breakEnt,hseq,bEnt)
#     aa     <- tmp$accept
#     if(aa == 1)accept['ent'] <- 0
#     if(aa == 0)accept['ent'] <- accept['ent'] + 1
#   }
#
#   if(pFlag){
#     tmp <- healthPars(bProp, priorBP, priorVBP, breakProp, breakPropLims, pStates, health, accept['prop'])
#     bProp      <- tmp$bgg
#     breakProp <- tmp$breakg
#     likeg <- likeg + tmp$loglik
#     #  plotLogit(breakPropLims, breakProp, hseq, bProp)
#     aa     <- tmp$accept
#     if(aa == 1)accept['prop'] <- 0
#     if(aa == 0)accept['prop'] <- accept['prop'] + 1
#   }

  if(rkFlag){
    tmp <- healthPars(bRake, priorBR, priorVBR, breakRake, breakRakeLims, rake, health, accept['rake'])
    bRake     <- tmp$bgg
    breakRake <- tmp$breakg
    likeg <- likeg + tmp$loglik
    # plotLogit(breakRakeLims,breakRake,hseq,bRake)
    aa     <- tmp$accept
    if(aa == 1)accept['rake'] <- 0
    if(aa == 0)accept['rake'] <- accept['rake'] + 1
  }

  if(skFlag){
    tmp <- healthPars(bSkin, priorBS, priorVBS, breakSkin, breakSkinLims, skin, health, accept['skin'])
    bSkin     <- tmp$bgg
    breakSkin <- tmp$breakg
    likeg <- likeg + tmp$loglik
    aa     <- tmp$accept
    if(aa == 1)accept['skin'] <- 0
    if(aa == 0)accept['skin'] <- accept['skin'] + 1
  }

  if(cyFlag){
    tmp <- healthPars(bCyam, priorBC, priorVBC, breakCyam, breakCyamLims, cyam, health, accept['cyam'])
    bCyam     <- tmp$bgg
    breakCyam <- tmp$breakg
    likeg <- likeg + tmp$loglik
    aa     <- tmp$accept
    if(aa == 1)accept['cyam'] <- 0
    if(aa == 0)accept['cyam'] <- accept['cyam'] + 1
  }

  if(caFlag){
    tmp <- healthPars(bCalf, priorBF, priorVBF, breakCalf, breakCalfLims, calves, health, accept['calf'])
    bCalf     <- tmp$bgg
    breakCalf <- tmp$breakg
    likeg <- likeg + tmp$loglik
    aa     <- tmp$accept
    if(aa == 1)accept['calf'] <- 0
    if(aa == 0)accept['calf'] <- accept['calf'] + 1
  }
  # }


  dev <- dev -2* likeg

  lamda <- updateLamda()

  aLam <- updateAlam()


  bggibbs[g,] <- bg
  bqgibbs[g,] <- bq
  if(bfFlag){bhgibbs[g,] <- bh}
  if(rkFlag){brgibbs[g,] <- bRake}
  if(skFlag){bsgibbs[g,] <- bSkin}
  if(cyFlag){bcgibbs[g,] <- bCyam}
  if(caFlag){bfgibbs[g,] <- bCalf}
#   if(entFlag){begibbs[g,] <- bEnt}
#   if(pFlag){bpgibbs[g,] <- bProp}

  hhgibbs[g,] <- breakBf
  lgibbs[g]  <- lamda[1]
  agibbs[g]  <- aLam
  suml       <- suml + lamda
  summ       <- summ + maleMoveProb
  sumf       <- sumf + femMoveProb
  sumx       <- sumx + xMoveProb
  summ2       <- summ2 + maleMoveProb^2
  sumf2       <- sumf2 + femMoveProb^2
  sumx2       <- sumx2 + xMoveProb^2

  sumLamda   <- sumLamda + lamda
  sumLamda2  <- sumLamda2 + lamda^2

  hh  <- health
  hhs <- matrix(death,n,nt)
  hh[htime > hhs] <- 0

  hh[is.na(hh)] <- 0
  sumh <- sumh + hh
  sumh2 <- sumh2 + hh^2

  # Getting the health of the population Sub-groups health tabulated
  # Population level health (Adult males and Juveniles)
  healthMale <- hh
  for(i in 1:nrow(m.idx)){
    healthMale[rownames(healthMale) == m.idx$SightingEGNo[i], 1:(m.idx$MinDateInt[i] - 1)] <- NA  # before first sight
  }
  healthMale[!rownames(healthMale) %in% m.idx$SightingEGNo,] <- NA
  newhgibbsMales  <- newhgibbsMales + healthMale
  newhgibbsMales2 <- newhgibbsMales2 + healthMale^2

  # 3.1 - Adult Males
  healthAdMale <- hh
  for(i in 1:nrow(am.idx)){
    healthAdMale[rownames(healthAdMale) == am.idx$SightingEGNo[i], 1:(am.idx$MinDateInt[i] - 1)] <- NA  # before first sight
  }
  healthAdMale[!rownames(healthAdMale) %in% am.idx$SightingEGNo,] <- NA
  newhgibbsAdMales  <- newhgibbsAdMales + healthAdMale
  newhgibbsAdMales2 <- newhgibbsAdMales2 + healthAdMale^2

 # 3.2 - Young Juveniles
  healthJuv <- hh
  for(i in 1:nrow(jd1.idx)){
    healthJuv[rownames(healthJuv) == jd1.idx$SightingEGNo[i], 1:(jd1.idx$MinDateInt[i] - 1)] <- NA
  }
  healthJuv[!rownames(healthJuv) %in% jd1.idx$SightingEGNo,] <- NA
  newhgibbsJuvs1 <- newhgibbsJuvs1 + healthJuv
  newhgibbsJuvs12 <- newhgibbsJuvs12 + healthJuv^2

 # 3.2 - Old Juveniles
  healthJuv <- hh
  for(i in 1:nrow(jd2.idx)){
    healthJuv[rownames(healthJuv) == jd2.idx$SightingEGNo[i], 1:(jd2.idx$MinDateInt[i] - 1)] <- NA
  }
  healthJuv[!rownames(healthJuv) %in% jd2.idx$SightingEGNo,] <- NA
  newhgibbsJuvs3 <- newhgibbsJuvs3 + healthJuv
  newhgibbsJuvs32 <- newhgibbsJuvs32 + healthJuv^2

# 4) Non-Reproductively Active Females (F > 9) # nullfem
  healthFem <- hh
  for(i in 1:nrow(nullfem)){
    healthFem[rownames(healthFem) == nullfem$EGNo[i], 1:(nullfem$MinDateInt[i] - 1)] <- NA
  }
  healthFem[!rownames(healthFem) %in% nullfem$EGNo,] <- NA
  newhgibbsNonRepFem  <- newhgibbsNonRepFem + healthFem
  newhgibbsNonRepFem2 <- newhgibbsNonRepFem2 + healthFem^2

  # 5) Delayed First Reproduction Females (F > 9) # delfem
  healthFem <- hh
  for(i in 1:nrow(delfem)){
    healthFem[rownames(healthFem) == delfem$EGNo[i], 1:(delfem$MinDateInt[i] - 1)] <- NA
  }
  healthFem[!rownames(healthFem) %in% delfem$EGNo,] <- NA
  newhgibbsDelRepFem  <- newhgibbsDelRepFem + healthFem
  newhgibbsDelRepFem2 <- newhgibbsDelRepFem2 + healthFem^2

  # 6) Pregnant females' health tabulated
  healthPregFem         <- hh
  healthPregFem[fp.idx] <- NA

  newhgibbsPregFem  <- newhgibbsPregFem  + healthPregFem
  newhgibbsPregFem2 <- newhgibbsPregFem2 + healthPregFem^2

  # 7) Lactating females' health tabulated
  healthLacFem         <- hh
  healthLacFem[fl.idx] <- NA

  newhgibbsLacFem  <- newhgibbsLacFem  + healthLacFem
  newhgibbsLacFem2 <- newhgibbsLacFem2 + healthLacFem^2

  # 8) Available to Calve females' health tabulated
  healthAvCavFem          <- hh
  healthAvCavFem[fac.idx] <- NA

  newhgibbsAvCavFem  <- newhgibbsAvCavFem  + healthAvCavFem
  newhgibbsAvCavFem2 <- newhgibbsAvCavFem2 + healthAvCavFem^2

  # 9) Resting females' health tabulated
  healthRestFem          <- hh
  healthRestFem[fr.idx]  <- NA

  newhgibbsRestFem  <- newhgibbsRestFem  + healthRestFem
  newhgibbsRestFem2 <- newhgibbsRestFem2 + healthRestFem^2


  if(g %in% c(50,100,200,500)){
    propBS     <- cov(bggibbs[1:(g-1),])
    bcovHealth <- cov(hhgibbs[1:(g-1),])
    diag(propBS)[diag(propBS) < .001] <- .001
  }

if(g==10000){save.image(file = paste(modelname,"_", g, '_wkspc.rdata', sep = ''))}

   if(RobHome){
    if(g==10000 & iter == 2){ngg <- g-1}
    if(g==10000 & iter == 2){save.image(file = paste(modelname,"_", g, '_wkspc.rdata', sep = ''))}
   }
}

close(pb)
