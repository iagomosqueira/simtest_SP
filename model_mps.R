# model_mps.R - DESC
# /home/mosquia/Active/Doing/SP_IOTC/simtest_SP/model_mps.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# - NO oem noise
#   - catch vs. vb

library(mse)

library(JABBA)
library(spict)
library(TMB)
library(SAMtool)

source("utilities.R")

# LOAD SWO om

load('bootstrap/data/swo.RData')

# BUG: MISSING in oem$idx
catch.wt(observations(oem)$idx[[1]]) <- stock.wt(stock(om))[, ac(1994:2046)]
catch.wt(observations(oem)$idx[[2]]) <- stock.wt(stock(om))[, ac(1994:2046)]

stock.wt(observations(oem)$stk)[,-1] <- stock.wt(observations(oem)$stk)[,1]

library(doParallel)
registerDoParallel(5)

# --- MP0

mseargs <- list(iy=2018, fy=2042, frq=3)

its <- which(unitSums(ssb(om)[, '2022']) == min(unitSums(ssb(om)[, '2022'])))
its <- seq(75)

om <- iter(om, its)
oem <- iter(oem, its)

# - JABBA

inits <- c(K=310000, r=0.42, q=c(0.000005263, 0.000005338))

control <- mpCtrl(list(
  est = mseCtrl(method=jabba.sa, args=list(inits=inits)),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$MSY) * 0.90,
      metric="Bdepletion", output="catch", dlow=0.85, dupp=1.15))))

control <- mpCtrl(list(
  est = mseCtrl(method=jabba.sa, args=list(inits=inits)),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$FMSY) * 0.90,
      min=0.001, metric="Bdepletion", output="fbar", dlow=0.85, dupp=1.15))))

system.time(
mp0jab <- mp(om, oem=oem, control=control, args=mseargs, parallel=TRUE)
)

# 48 min @bruce w/ 50 cores

# - DLMSP

# MSY

control <- mpCtrl(list(
  est = mseCtrl(method=dlmsp.sa, args=list(n_seas=1L, random="log_B_dev",
    prior_dist=list(r=c(log(0.42), 0.4), MSY=c(log(30000), 0.3)))),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$MSY) * 0.40,
      metric="Bdepletion", output="catch", dlow=0.85, dupp=1.15))))

# FMSY

control <- mpCtrl(list(
  est = mseCtrl(method=dlmsp.sa, args=list(n_seas=1L, random="log_B_dev",
    prior_dist=list(r=c(log(0.42), 0.4), MSY=c(log(30000), 0.3)))),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$FMSY) * 0.90,
      metric="Bdepletion", output="fbar"))))

# 12 min @bruce w/ 50 cores

system.time(
mp0dlm <- mp(om, oem=oem, control=control, args=mseargs, parallel=TRUE, 
  verbose=TRUE)
)

plot(FLStocks(OM=nounit(window(stock(om), end=2018)),
  MP=nounit(stock(mp0dlm))))

# tracking
tracking(mp0dlm)

# CHECK model convergence
iterMeans(tracking(mp0dlm)['conv.est',])


# TUNE

mets <- list(F=function(x) unitMeans(fbar(x)),
  SB=function(x) unitSums(ssb(x)))

mp0tun <- tunebisect(om, oem, control=control,
  tune=list(target=mean(refpts(om)$FMSY) * c(0.10, 10)), prob=0.6,
  metrics=mets, statistic=statistics['green'], 
  args=mseargs, years=seq(2034, 2038))


# COMPUTE performance
per <- performance(mp0tun, statistics['green'], years=list(2034:2038),
  metrics=mets)

# CHECK tuning performance
per[year == 2038 & statistic=='green', .(pgreen=mean(data)), by=mp]


plot(FLStocks(OM=nounit(window(stock(om), end=2018)),
  MP1=nounit(stock(mp0tun[[1]])),
  MP2=nounit(stock(mp0tun[[2]]))
  ))






# RUN for om iter w/min(ssb).


# - SPICT

control <- mpCtrl(list(
  est = mseCtrl(method=spict.sa),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$MSY) * 0.90,
      metric="Bdepletion", output="catch", dlow=0.85, dupp=1.15))))

system.time(
mp0spi <- mp(om, oem=oem, control=control, args=mseargs, parallel=FALSE, 
  verbose=TRUE)
)


plot(om, JABBA=mp0jab, DLMSP=mp0dlm, SPICT=mp0spi)

plot(FLStocks(OM=nounit(window(stock(om), end=2018)),
  MP=nounit(stock(mp0jab))))
