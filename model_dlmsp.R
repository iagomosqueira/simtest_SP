# model_dlmsp.R - DESC
# test_SP/model_dlmsp.R

# Copyright (c) WUR, FAO 2023.
# Author: Henning WINKER (FAO) <iago.mosqueira@wur.nl>
#         Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(TMB)
library(SAMtool)
load_all('~/FLR/pkgs/MSE_PKG@flr/mse')
library(mse)

source("utilities.R")

# - SIM {{{

# LOAD OM

load("data/sims.RData")

# SET inputs

args <- list(it=dim(om)[6], ay=2020)

tracking <- FLQuant(dimnames=list(metric="conv.est", year=1951:2020,
  iter=seq(args$it)))

# R01: No priors

r01 <- dlmsp.sa(om, idx=FLIndices(A=idx), args, tracking)$ind

# R02: r and MSY priors

r02 <- dlmsp.sa(om, idx=FLIndices(A=idx), args, tracking,
  prior_dist=list(r=c(0.4, 0.2), MSY=c(log(50), 0.4)))$ind

# R03: r and MSY priors, 4 seasons

r03 <- dlmsp.sa(om, idx=FLIndices(A=idx), args, tracking,
  prior_dist=list(r=c(0.4, 0.2), MSY=c(log(50), 0.4)), n_seas=4L)$ind

# r04: random "log_B_dev"

r04 <- dlmsp.sa(om, idx=FLIndices(A=idx), args, tracking,
  prior_dist=list(r=c(0.4, 0.2), MSY=c(log(50), 0.4)), n_seas=4L,
  random="log_B_dev")$ind

# r05: no state_space

r05 <- dlmsp.sa(om, idx=FLIndices(A=idx), args, tracking,
  prior_dist=list(r=c(0.4, 0.2), MSY=c(log(50), 0.4)), n_seas=4L,
  random="log_B_dev", state_space=FALSE)$ind

# OM

rom <- FLQuants(F=fbar(om), Fstatus=fbar(om) / refpts(om)$Fmsy, B=tsb(om),
  Bstatus=tsb(om) / refpts(om)$Btgt, Bdepletion=tsb(om) / refpts(om)$B0)

# - RESULTS

rs <- list(OM=rom, R01=r01, R02=r02, R03=r03, R04=r04, R05=r05)

res <- lapply(setNames(nm=names(r01)), function(x)
  FLQuants(lapply(rs, "[[", x)))

pls <- Map(function(x,y) plot(x) + ggtitle(y) +ylim(c(0,NA)),
  x=res, y=names(res))

Reduce('+', pls)

save(res, pls, file="model/dlmsp.RData", compress="xz")

# }}}

# - SWO IOTC MP {{{

load('~/Active/Old/IOTC/SWO_MSE@iotc/swo/OM/output/om.Rdata')

library(doParallel)
registerDoParallel(2)

mseargs <- list(iy=2018, fy=2035, frq=3)

om <- iter(om, seq(100))
oem <- iter(oem, seq(100))

# MP0

control <- mpCtrl(list(
  est = mseCtrl(method=dlmsp.sa, args=list(n_seas=1L, random="log_B_dev",
    prior_dist=list(r=c(log(0.3, 0.2)), MSY=c(log(30000), 0.3)))),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$MSY) * 0.90,
      metric="Bdepletion", output="catch", dlow=0.85, dupp=1.15))))

system.time(
mp0 <- mp(om, oem=oem, control=control, args=mseargs, parallel=TRUE)
)

# MP1
args(control$est)$random <- NULL

system.time(
mp1 <- mp(om, oem=oem, control=control, args=mseargs, parallel=TRUE)
)

# MP2
args(control$est)$n_seas <- 4L

system.time(
mp2 <- mp(om, oem=oem, control=control, args=mseargs, parallel=TRUE)
)


save(mp0, mp1, mp2, file="model/swo_dlsmp.RData", compress="xz")


# PLOTS

# PLOT OM + MP run
plot(FLStocks(OM=nounit(window(stock(om), end=2018)),
  MP=nounit(stock(mp0))))

plot(tracking(mp0)[c("B.om")] / refpts(om)$B0,
  tracking(mp0)[c("met.hcr")])


# }}}

# QUESTIONS

# - state space
# - B, SB, VB
# - start
# - I_lambda

# Schaefer, r prior lnrom prior 0.2-0.3, CV 0.5 and process error (if fixed) around 0.07-0.1


# FIX prior, variance. RUN dlmsp, apply with Jabba
# SCALE K prior

