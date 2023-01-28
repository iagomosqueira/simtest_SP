# data.R - Builds an OM used for testing SP models.
# test_SP/data.R

# Copyright (c) WUR, FAO 2023.
# Author: Henning WINKER (FAO) <iago.mosqueira@wur.nl>
#         Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

library(FishLife)
library(SPMpriors)
library(FLRef)
library(FLife)
library(JABBA)
library(patchwork)


# --- Life history

# A set of parameters are sampled from FishLife

alb <- flmvn_traits(Genus="Thunnus", Species="alalunga",
  h=c(0.6,0.9), Plot=FALSE)$traits

# Other parameters are specified

lh <- list(linf = alb[alb$trait=="Loo","mu.stk"],
  k =  alb[alb$trait=="K","mu.stk"],
  t0 = -0.5,
  tm = alb[alb$trait=="tm","mu.stk"],
  a = 0.00001,
  b = 3.04,
  tmax =  ceiling(alb[alb$trait=="tmax","mu.stk"]),
  s = alb[alb$trait=="h","mu.stk"])

# Check the corresponding Von Bertalannfy growth curve

lena <- with(lh, vonbert(linf, k, t0, age=0:tmax))

# M using Then (), with tmax1 = 22 and tmax2 = 12

lh$M <- 4.899 * lh$tmax ^- 0.916

#  Leslie matrix `r`

bio <- with(lh, jbio(amax=tmax, nsexes=1, Loo=linf, k=k, t0=t0, aW=1000*a, bW=b,
  mat=c(tm, tm*1.2, 1), M=M, h=s))

r <- jbleslie(bio)$r


# --- Generic OM for ALB

# Life-history parameter set com,pleted by lhPar

par <- with(lh, lhPar(FLPar(linf=linf, k=k, t0=t0, a=a, b=b, a50=2, s=s, m1=M)))

# Constructor a basic `FLBRP()` (equilibrium) object.

eql <- lhEql(par, range = c(min = 0, max = lh$tmax,
  minfbar = 1, maxfbar = lh$tmax-1, plusgroup=lh$tmax))

# Adjusts the inbuilt Gislason M to a scaled Lorenzen $M$

m(eql)[] <- lorenzen(stock.wt(eql), Mref=lh$M, Aref=2)

 # Specifying a new maturity ogive

mat(eql) <- newselex(mat(eql), FLPar(S50=lh$tm, mat95=lh$tm * 1.5, Smax=1000,
  Dcv=0.5,Dmin=0.1))

# Updates the stock recruitment steepness if biological parameters are modified

params(eql)

eql <- updsr(eql, s=lh$s) 

params(eql)

# Time horizon can be adjusted with `fbar` range, 1981-2020

f0 <- 0.01

fbar(eql) <- FLQuant(f0, dimnames=list(year=1951:2020))

# Convert FLBRP into FLStock + SRR

stk <- as(eql,"FLStock")
units(stk) <- standardUnits(stk)

srr <- as(as(eql,"predictModel"),"FLSR")

# Set new selectivity parameters

selpar <- FLPar(
  # catch before maturity
  S50=3, S95=4.2,
  # Domed-shaped
  Smax=7, Dcv=0.5, Dmin=0.2)

# Modify selectivity

harvest(stk) <-  f0 * newselex(harvest(stk), selpar)

# Plots weight-at-age, maturity-at-age, natural mortality and selective

ggplot(FLQuants(stk,"m","catch.sel","mat","catch.wt"))+
  geom_line(aes(age,data))+
  facet_wrap(~qname,scale="free")+theme_bw()


# --- Estimating priors from FLStock and SRR

#' First, it is straight forward to compute the intrinsic rate of population increase from a Leslie for a specified steepness value. However, this $r$ estimate ignores selectivity and should only be used in the context of a Schaefer model with $F_{MSY} = r/2$.  

r.leslie <- mean(productivity(stk, s=lh$s)$r)

# compare
c(r=r, r.leslie=r.leslie, FMSY=r.leslie / 2)

#' An alternative is to estimate $r$ and the shape parameter $n$ as function $MSY$, $VB_{MSY}$ and $VB_0$ from an age-structured equilibrium model (ASEM), where $BB$ denotes that the vulnarable (or exploitable) biomass as function of selectivity. 

brp <- brp(FLBRP(stk, srr))

r.pella <- asem2spm(brp)
r.pella

plotpf(brp, rel=FALSE) + plotpf(brp, rel=TRUE) + plot_layout(guides = "collect")


# --- Simulating stock dynamics with evolutionary F-trajectories

# Estimate refpts
brp <- computeFbrp(stk, srr, proxy="msy", blim=0.3, type="btgt")
stk <- FLStockR(stk)
refpts(stk) <- Fbrp(brp)

ploteq(brp)

its <- 100 

# SSB0 for brp
b0 <- an(refpts(eql)["virgin","ssb"])

# one-way downhill projection
control = FLPar(Feq=0.15,Frate=0.05,Fsigma=0.15,SB0=b0,
  minyear=dims(stk)$minyear+1,
  maxyear=dims(stk)$maxyear,its=its)

# Forecasted with random recruitment under an evolving F-trajectory

set.seed(1234)
# Random recruitment deviations with sigR = 0.5 and AR1 rho = 0.3
rec_devs = ar1rlnorm(0.3, 1951:2020, its, 0, 0.5)
# propagate desired iterations
stki <- propagate(stk, its)
# create OM 
om <- rffwd(stki, srr, control=control, deviances=rec_devs)
refpts(om) <- refpts(stk)


plotAdvice(om)

# Generate LL CPUE, flat-toped selectivity

cpue.sel <- newselex(harvest(stk),
  FLPar(S50=3.5, mat95=4.2, Smax=100, Dcv=0.5, Dmin=0.05))

idx <- window(bioidx.sim(om,sel=cpue.sel,sigma=0.25,q=0.001),start=1970)

flqs= FLQuants(
  "B/Bmsy"=ssb(om)/om@refpts["Btgt"],
  "F/Fmsy"=fbar(om)/om@refpts["Fmsy"],
  "Catch/MSY"=catch(om)/om@refpts["Yeq"],
  Index=idx@index%/%yearMeans(idx@index))

ggplot(iter(flqs,1:20))+
  geom_line(aes(year,data,col=ac(iter)),alpha=0.5)+
  theme_bw()+facet_wrap(~qname)+
  theme(legend.position = "none")+geom_hline(yintercept = 1,linetype=2)+
  scale_color_manual(values=c(rainbow(20)))

#

save(om, idx, file="data/om.RData", compress="xz")
