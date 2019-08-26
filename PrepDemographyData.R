library(ggplot2)
library(lubridate)
library(stringr)
library(gdata)
library(dplyr)
# Demography Data Prep
rm(list=ls())
source(file='/Users/rob/Documents/code/rss10/rightwhales/SPUE.Sim.R')
source(file='/Users/rob/Documents/code/rss10/rightwhales/MoveProb.R')
# source(file='/Users/rob/Documents/code/rss10/rightwhales/makeTangle.r')
wkdir <- "/Users/rob/Documents/research/projects/PCAD/rightwhales"
datadir <- '/data'
setwd(wkdir)
load(file = "data/egsightings_2016-04-29.rdata")
name <- paste("egDemogData_", Sys.Date(), ".rdata", sep = '') # output file name

batch$PhotoCount <- batch$AerialCount + batch$ShipboardCount
batch$Date <- ymd(paste(batch$EndYear, batch$EndMonth, 1, sep = '/'))
# sights <- sights[-which(sights[,'SightingMonth']==0),]
sights[which(sights[, 'SightingDay'] == 0), 'SightingDay'] <- 1
sights$Date       <- ymd(paste(sights$SightingYear, sights$SightingMonth, sights$SightingDay, sep = '/'))
calves            <- read.csv(file = "data/CalvingEventsTable_2014-05-07.csv")

# Assemble a vector to concatenate month-years
first.date       <- min(sights$Date,na.rm=T)
last.date        <- max(sights$Date,na.rm=T)
moyr             <- seq.Date(as.Date(first.date), as.Date(last.date), by = '1 months')
sights$monthyear <- findInterval(as.Date(sights$Date), moyr)

# Get the longitude on the right side of the meridian
sights$LongGood  <- -1*sights$Longitude

# Make a summary of the sights:
# sights.sum: sightings per animal per region
# calves.sum: # calves per animal
sights.sum <- sights %>%
  group_by(SightingEGNo, RegionCode, monthyear) %>%
  summarise(nsights = n())

calves.sum <- calves %>%
  group_by(EGNo) %>%
  summarise(num.calves = n())


effort <- spue.sim(moyr , yr = unique(sights$SightingYear))
nextmat <- MoveProb(moyr , yr = unique(sights$SightingYear))

# I'm doing this to get the narrative movement transition values from Philip Hamilton at NEAq incorporated
# Sep 14, 2016, adding in Philip's new and revised priors
adultMales    <- read.csv(file = 'data/adultMales_PhilipOriginal_Nov2011.csv', header = T)
adultMalesNew <- read.csv(file = 'data/adultMales_PhilipUpdated_Sep2016.csv', header = T)

adultFemales    <- read.csv(file = 'data/adultFemales_PhilipOriginal_Nov2011.csv', header = T)
adultFemalesNew <- read.csv(file = 'data/adultFemales_PhilipUpdated_Sep2016.csv', header = T)

xGender    <- read.csv(file = 'data/xGender_PhilipOriginal_Nov2011.csv', header = T)
xGenderNew <- read.csv(file = 'data/xGender_PhilipUpdated_Sep2016.csv', header = T)

juvMales    <- read.csv(file = 'data/juvMales_PhilipOriginal_Nov2011.csv', header = T)
juvMalesNew <- read.csv(file = 'data/juvMales_PhilipUpdated_Sep2016.csv', header = T)

juvFemales    <- read.csv(file = 'data/juvFemales_PhilipOriginal_Nov2011.csv', header = T)
juvFemalesNew <- read.csv(file = 'data/juvFemales_PhilipUpdated_Sep2016.csv', header = T)


# Save the data into one rdata file
save(batch, calves.sum, sights.sum, sights, calves, effort, nextmat, adultMales, adultFemales, xGender, juvMales, juvFemales, file = name)
