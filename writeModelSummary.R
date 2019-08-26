zz <- file(paste('gibbsoutput/', modelname, "modelParameterSum.txt", sep = ''), "w")  
cat('Model name = ', modelname , file = zz,  "\n")
cat('# Gibbs Steps = ', ng , file = zz,  "\n")

cat("True/False Flags", file = zz, '\n')
cat('subAn = ', subAn,  file = zz,  "\n")
cat('subYrs = ', subYrs,  file = zz,  "\n")
cat('FixSurv = ', FixSurv,  file = zz,  "\n")
cat('Survhack = ', Survhack,  file = zz,  "\n")
cat('Use Regional Intercept (REGINTERCEPT) = ', REGINTERCEPT,  file = zz,  "\n")
cat('Impute missing data (missDataImp) = ', missDataImp,  file = zz,  "\n")
cat('Movement Prior Sensitivity Analyis (movePsa) = ', movePsa, file = zz, '\n')
if(movePsa){cat('Movement Prior Sensitivity Value (movePsaVal) = ', movePsaVal, file = zz, '\n')}

cat("Parameter Values\n", file = zz)
cat('Visual Health Parameters = ', hnamevec , file = zz,  "\n")
cat('Amount proposed health can change each month (hoffset)', hoffset, file = zz, '\n')
if(bfFlag){cat('Body Fat Limits = ', breakBfLims , file = zz,  "\n")}
if(rkFlag){cat('Rake Mark Limits = ', breakRakeLims , file = zz,  "\n")}
if(skFlag){cat('Skin Limits = ', breakSkinLims , file = zz,  "\n")}
if(cyFlag){cat('Cyamid Limits = ', breakCyamLims , file = zz,  "\n")}
cat('Health Process Variance (muh) = ', muh , file = zz,  "\n")

close(zz)
