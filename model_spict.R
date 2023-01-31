# model_spict.R - DESC
# /home/mosquia/Active/Doing/SP_IOTC/test_SP/model_spict.R

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

load("data/om.RData")

# SET dlmp.sa inputs

args <- list(it=dim(om)[6], ay=2020)

tracking <- FLQuant(dimnames=list(metric="conv.est", year=1951:2020,
  iter=seq(args$it)))

# CALL dlmp.sa w/priors

system.time(
s01 <- spict.sa(om, idx=FLIndices(A=idx), args, tracking)
)

# PLOTS

# estimates by iter
ggplot(s01$ind, aes(x=year,y=data,group=iter)) + ylim(c(0, NA)) +
  geom_line() + facet_wrap(~qname, scales="free")

# B ~ B(om)
plot(s01$ind$B, stock(om), ssb(om)) + ggtitle("B, SSB, Bhat")

# B ~ VB(om)
plot(s01$ind$B, vb(window(om, start=1970), sel.pattern(idx))) +
  ggtitle("VB, Bhat")

# }}}


# --- SWO

load("bootstrap/data/swo.RData")

stk <- iter(nounit(swo), seq(500))
idx <- iter(swoidx, seq(500))

args <- list(it=dims(stk)$iter, ay=dims(stk)$endyear)
tracking <- FLQuant(dimnames=list(metric="conv.est", year=1951:2020,
  iter=seq(args$it)))

system.time(
sp01 <- spict.sa(stk, idx, args=args, tracking=tracking)
)

plot(sp01$ind)

ggplot(sp01$ind, aes(x=year, y=data, group=iter)) +
  geom_line() +
  facet_grid(qname~., scales="free")

ggplot(FLQuants(SP=sp01$ind$B, OM=stock(swo))) +
  geom_line(aes(x=year, y=data, group=iter)) +
  facet_grid(qname~., scales="free")

