# Progress Bar
pb <- txtProgressBar(min = 0, max = ng, style = 3)
# dev <- 0

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
