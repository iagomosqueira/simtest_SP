# model_dlmsp.R - DESC
# /home/mosquia/Active/Doing/SP_IOTC/test_SP/model_dlmsp.R

# Copyright (c) WUR, FAO 2023.
# Author: Henning WINKER (FAO) <iago.mosqueira@wur.nl>
#         Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(TMB)
library(SAMtool)
library(mse)

source("utilities.R")


# - SIM {{{

# LOAD OM

load("data/om.RData")

# SET dlmp.sa inputs

args <- list(it=dim(om)[6])

tracking <- FLQuant(dimnames=list(metric="conv.est", year=1951:2020,
  iter=seq(args$it)))

# CALL dlmp.sa w/priors

s01 <- dlmsp.sa(om, idx=FLIndices(A=idx), args, tracking,
  prior_dist=list(r=c(0.4, 0.2), MSY=c(log(50), 0.4)))

# PLOTS

# estimates by iter
ggplot(s01$ind, aes(x=year,y=data,group=iter)) + ylim(c(0, NA)) +
  geom_line() + facet_wrap(~qname, scales="free")

# B ~ B(om)
plot(s01$ind$B, stock(om))

# B ~ VB(om)
plot(s01$ind$B, vb(window(om, start=1970), sel.pattern(idx)))

# }}}


# - SWO IOTC {{{

load("bootstrap/data/swo.RData")

args <- list(it=dim(swo)[6])

tracking <- FLQuant(dimnames=list(quant="conv", year=2018, iter=seq(args$it)))

# - SWO OM

tes <- dlmsp.sa(swo, swoidx, args, tracking)




plot(tes$ind) + ylim(c(0, NA))

ggplot(tes$ind, aes(x=year,y=data,group=iter)) + ylim(c(0, NA)) +
  geom_line() + facet_wrap(~qname, scales="free")

ggplot(FLQuants(C=catch(stk), I1=index(idx[[1]]), I2=index(idx[[2]])),
  aes(x=year,y=data, group=iter, colour=qname)) +
  ylim(c(0, NA)) +
  geom_line() + facet_wrap(~qname, scales="free")

ggplot(indices, aes(x=year,y=data,group=iter)) + ylim(c(0, NA)) +
  geom_line() + facet_wrap(~qname, scales="free")

# TODO: CHECK indices

# }}}

# QUESTIONS

# - state space
# - B, SB, VB
# - start
# - I_lambda

# Schaefer, r prior lnrom prior 0.2-0.3, CV 0.5 and process error (if fixed) around 0.07-0.1


# FIX prior, variance. RUN dlmsp, apply with Jabba
# SCALE K prior

