# model_jabba.R - DESC
# /home/mosquia/Active/Doing/SP_IOTC/test_SP/model_jabba.R

# Copyright (c) WUR, FAO 2023.
# Author: Henning WINKER (FAO) <iago.mosqueira@wur.nl>
#         Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(JABBA)
library(mse)

source("utilities.R")


# - SIM {{{

# LOAD OM

load("data/sims.RData")

# SET inputs

args <- list(it=dim(om)[6], ay=2020, y0=1951)

tracking <- FLQuant(dimnames=list(metric="conv.est", year=1951:2020,
  iter=seq(args$it)))

system.time(
r01 <- jabba.sa(om, idx=FLIndices(A=idx), args, tracking)[["ind"]]
)

# 100 iters, min

# OM

rom <- FLQuants(B=tsb(om), Bstatus=tsb(om) / refpts(om)$Btgt,
  Bdepletion=tsb(om) / refpts(om)$B0)

# - RESULTS

rs <- list(OM=rom, JABBA=r01)

res <- lapply(setNames(names(rom), nm=c("B", "B/B[MSY]", "B/B[0]")),
  function(x) FLQuants(lapply(rs, "[[", x)))

pls <- Map(function(x,y) plot(x) + ggtitle(parse(text=y)) + ylim(c(0,NA)),
  x=res, y=names(res))

Reduce('+', pls)

pubpng("report/sim_b.png", width=2200,
  Reduce('+', pls))

save(res, pls, file="model/jabba.RData", compress="xz")

# }}}

# - SWO IOTC MP {{{

load('~/Active/SWO_MSE@iotc/swo/OM/output/om.Rdata')

library(doParallel)
registerDoParallel(2)

args <- list(it=1, ay=2020, y0=1950)

tracking <- FLQuant(dimnames=list(metric="met", year=1951:2020,
  iter=seq(args$it)))

stk <- iter(window(nounit(stock(om)), end=2018), 1)

run <- jabba.sa(stk,
  idx=iter(window(observations(oem)$idx, end=2018), 1), args, tracking)

plot(stock(stk) / iter(refpts(om)$B0, 1))
plot(stock(stk) / iter(refpts(om)$B0, 1), run$ind$Bdepletion)

# MP0

mseargs <- list(iy=2018, fy=2035, frq=3)

om <- iter(om, seq(2))
oem <- iter(oem, seq(2))


control <- mpCtrl(list(
  est = mseCtrl(method=jabba.sa),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$MSY) * 0.90,
      metric="Bdepletion", output="catch", dlow=0.85, dupp=1.15))))

system.time(
mp0 <- mp(om, oem=oem, control=control, args=mseargs, parallel=TRUE)
)

tracking(mp0)['pid',]


# MP1
args(control$est)$random <- NULL

system.time(
mp1 <- mp(om, oem=oem, control=control, args=mseargs, parallel=5)
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

