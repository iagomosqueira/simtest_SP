# model_spict.R - DESC
# test_SP/model_spict.R

# Copyright (c) WUR, FAO 2023.
# Author: Henning WINKER (FAO) <iago.mosqueira@wur.nl>
#         Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(spict)
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

r01 <- spict.sa(om, idx=FLIndices(A=idx), args, tracking)

# OM

rom <- FLQuants(F=fbar(om), Fstatus=fbar(om) / refpts(om)$Fmsy, B=tsb(om),
  Bstatus=tsb(om) / refpts(om)$Btgt, Bdepletion=tsb(om) / refpts(om)$B0)

# - RESULTS

rs <- list(OM=rom, R01=r01$ind)

res <- lapply(setNames(nm=names(rom)), function(x)
  FLQuants(lapply(rs, "[[", x)))

pls <- Map(function(x,y) plot(x) + ggtitle(y) +ylim(c(0,NA)),
  x=res, y=names(res))

Reduce('+', pls)

save(res, pls, file="model/dlmsp.RData", compress="xz")

# }}}

# - SWO IOTC MP {{{

load('bootstrap/data/swo.RData')

mseargs <- list(iy=2018, fy=2035, frq=3)

om <- iter(om, seq(100))
oem <- iter(oem, seq(100))

# MP0

control <- mpCtrl(list(
  est = mseCtrl(method=spict.sa),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$MSY) * 0.90,
      metric="Bdepletion", output="catch", dlow=0.85, dupp=1.15))))

# BUG: iter 93, lowest biomass
args(projection(om)) <- list(maxF=0.5)

system.time(
mp0 <- mp(om, oem=oem, control=control, args=mseargs, parallel=FALSE)
)

# - PLOTS

# PLOT OM + MP run

plot(FLStocks(OM=nounit(window(stock(om), end=2018)),
  MP=nounit(stock(mp0))))

plot(tracking(mp0)[c("B.om")] / refpts(om)$B0,
  tracking(mp0)[c("met.hcr")])

# }}}
