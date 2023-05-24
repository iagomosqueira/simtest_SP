# model_swo.R - DESC
# /home/mosquia/Active/Doing/SP_IOTC/simtest_SP/model_swo.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)
library(JABBA)
library(spict)
library(TMB)
library(SAMtool)

source("utilities.R")

# LOAD SWO om

load('bootstrap/data/swo.RData')


# --- COMPARE SPs on SWO stk 2018

its <- seq(250)
its <- sample(seq(500), 2)
its <- c(1, 10)
its <- 500

args <- list(it=length(its), ay=2019, data_lag=1, y0=1950)

tracking <- FLQuant(dimnames=list(metric="met", year=1951:2046,
  iter=its))

stk <- iter(window(nounit(stock(om)), end=2018), its)
idx <- iter(window(observations(oem)$idx, end=2018), its)

# - JABBA

# TODO: USE 35 quantile depletion Jabba

inits <- c(K=310000, r=0.42, q=c(0.000005263, 0.000005338))

rjab <- jabba.sa(stk, idx, args, tracking, model.type="Schaefer",
  inits=inits)

nrjab <- jabba.sa(stk, idx, args, tracking, model.type="Schaefer")

# - SPICT

rspic <- spict.sa(stk, idx=idx, args, tracking)

# - DLMSP

rdlm <- dlmsp.sa(stk, idx=idx, n_seas=1L, random="log_B_dev",
  prior_dist=list(r=c(log(0.42), 0.4), MSY=c(log(30000), 0.3)),
  args, tracking)

save(stk, idx, rjab, rspic, rdlm, file='model_swo_sa.RData')


# ---

# - PROJECT om

ctrl <-  fwdControl(lapply(2023:2035, function(x)
  list(year=x, quant='fbar', value=refpts(om)$FMSY * 1.25)))

om <- fwd(om, control=ctrl)

# - GENERATE indices for new years

catch.wt(observations(oem)$idx[[1]]) <- stock.wt(stock(om))[, ac(1994:2046)]
catch.wt(observations(oem)$idx[[2]]) <- stock.wt(stock(om))[, ac(1994:2046)]

observations(oem)$idx[[1]][, ac(2023:2035)] <- survey(
  stock(om)[, ac(2023:2035)], observations(oem)$idx[[1]][, ac(2023:2035)])
observations(oem)$idx[[2]][, ac(2023:2035)] <- survey(
  stock(om)[, ac(2023:2035)], observations(oem)$idx[[2]][, ac(2023:2035)])

idx <- window(observations(oem)$idx, end=2035)

args <- list(it=500, ay=2019, data_lag=1, y0=1950)

tracking <- FLQuant(dimnames=list(metric="met", year=1951:2046,
  iter=seq(500)))

inits <- c(K=310000, r=0.42, q=c(0.000005263, 0.000005338))

frjab <- jabba.sa(window(nounit(stock(om)), end=2035), idx,
  args, tracking, model.type="Schaefer", inits=inits)







# PDF

stk <- window(nounit(stock(om)), end=2018)

pdf("report/swo_fits.pdf")

plot(FLQuants(VB=vb(stk, sel=propagate(sel.pattern(idx[[1]])[,1,,,,1], 500)),
  JABBA=rjab$ind$B,
  SPICT=rspic$ind$B,
  DLMSP=rdlm$ind$B)) + ggtitle("Vulnerable biomass") +
  ylim(c(0, NA))

plot(FLQuants(B=stock(stk),
  JABBA=rjab$ind$B,
  SPICT=rspic$ind$B,
  DLMSP=rdlm$ind$B)) + ggtitle("Total biomass") +
  ylim(c(0, NA))

plot(FLQuants(OM=stock(stk) %/% refpts(om)['B0',],
  JABBA=rjab$ind$Bdepletion,
  SPICT=rspic$ind$Bdepletion,
  DLMSP=rdlm$ind$Bdepletion)) + ggtitle("Depletion") +
  ylim(c(0,NA))

dephat <- FLQuants(OM=stock(stk)[, ac(1950:2018)] %/%
    refpts(om)['B0',],
  JABBA=rjab$ind$Bdepletion,
  SPICT=rspic$ind$Bdepletion, 
  DLMSP=rdlm$ind$Bdepletion)

dat <- data.table(as.data.frame(window(dephat, start=2018, end=2018))[, c(7,8)])

library(ggpmisc)

plot(dephat) +
  annotate("table", x=1960, y=0.25,
  label=dat[, .(mean=mean(data), SD=sd(data)), by=qname]) +
  ylim(c(0,NA))

ggplot(dephat, aes(x=year, y=data, group=iter)) +
  geom_line() + ylim(c(0,NA))+ ggtitle("Depletion") +
  ylab("") + xlab("") +
  annotate("table", x=1960, y=0.25,
  label=dat[, .(mean=mean(data), SD=sd(data)), by=qname]) +
  facet_grid(qname~.) +
  ylim(c(0,NA))

dev.off()

