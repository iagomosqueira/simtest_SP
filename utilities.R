# utilities.R - DESC
# /home/mosquia/Active/Doing/SP_IOTC/test_SP/utilities.R

# Copyright (c) WUR, FAO 2023.
# Author: Henning WINKER (FAO) <iago.mosqueira@wur.nl>
#         Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# ---- dlmsp

# buildTMBdataDLMSP {{{

buildTMBdataDLMSP <- function(catch, indices, rescale=1, fix_sigma=1,
  state_space=TRUE, early_dev=c("all", "index"),
  prior_dist=list(r=c(NA, NA), MSY=c(NA, NA)), n_seas=4L) {

  # PARSE and CHECK arguments
  
  ny <- dim(catch)[2]

  # - COMPLETE catch series

  # - random

  # EXTEND indices to match catch
  indices <- window(indices, start=dims(catch)$minyear,
    end=dims(catch)$maxyear)

  # PARSE data list elements

  # C_hist
  C_hist <- c(catch)

  # I_hist, I_sd, I_lambda
  I_hist <- unname(do.call('cbind', lapply(indices, '[', drop=TRUE)))
  I_sd <- I_hist
  # TODO: PARSE I_sd
  I_sd[] <- 0.3
  nsurvey <- ncol(I_hist)
  I_lambda <- rep(1, nsurvey)

  n_itF <- 3L
  
  # PARSE state_space & early_dev, changes est_B_dev
  if(state_space) {
    early_dev <- match.arg(early_dev)
    if(early_dev == "all")
      est_B_dev <- rep(1, ny)
    if(early_dev == "index")
      est_B_dev <- rep(1, ny)
  } else {
    # 0 until fist index value (caross indices), then 1
    ipos <- min(apply(!is.na(I_hist), 2, which.max))
    est_B_dev <- c(rep(0, ipos), rep(1, ny - ipos))
  }

  # BUILD priors arguments
  use_prior <- !is.na(prior_dist)

  if(any(use_prior)) {
    prior_dist <- matrix(unlist(prior_dist), ncol=2, nrow=2,
      byrow=TRUE, dimnames=list(c("r", "MSY"), c("par1", "par2")))
  }

  # BUILD data list

  data <- list(model = "SP", C_hist = C_hist, rescale = rescale,
    I_hist = I_hist, I_sd = I_sd, I_lambda = I_lambda,
    fix_sigma = as.integer(fix_sigma),
    nsurvey = nsurvey, ny = length(C_hist), nstep = n_seas, 
    dt = 1/n_seas, n_itF = n_itF, est_B_dev = est_B_dev,
    sim_process_error = 0L,
    use_prior = use_prior, prior_dist = prior_dist)

  return(data)
}
# }}}

# dlmsp.sa {{{

dlmsp.sa <- function(stk, idx, args, tracking,
  start=list(), state_space=TRUE, early_dev=c("all", "index"),
  prior_dist=list(r=c(NA, NA), MSY=c(NA, NA)),
  n_seas=4L, rescale=1, random=NULL, verbose=FALSE) {
 
  # EXTRACT inputs
  catch <- catch(stk)
  indices <- lapply(idx, index)

  # DIMS
  ny <- dim(catch)[2]
  nsurvey <- length(indices)

  # TODO: CHECK inputs

  # OUTPUT
  empty <- catch %=% 0
  out <- vector(mode="list", length=args$it)

  # SETUP map
  map <- list(log_dep = factor(NA), log_n = factor(NA),
    log_tau = factor(NA), log_sigma = factor(rep(NA, nsurvey)), 
    log_B_dev = factor(rep(NA, ny)))

  # LOOP over iters
  for (i in seq(args$it)) {

    print(i)

    # EXTRACT single iter
    cab <- iter(catch, i)
    ind <- iter(indices, i)

    data <- buildTMBdataDLMSP(cab, ind, rescale=rescale,
      state_space=state_space, early_dev=early_dev,
      prior_dist=prior_dist, n_seas=n_seas)

    # SETUP params
    params <- list(
      # log_FMSY
      log_FMSY = log(0.20),
      # MSYx
      MSYx = log(mean(3 * c(cab) * rescale)),
      # log_dep
      log_dep = log(1),
      # log_n
      log_n = log(2),
      # log_sigma
      log_sigma = rep(log(0.05), nsurvey),
      # log_tau
      log_tau = log(0.1),
      # log_B_dev
      log_B_dev = rep(0, ny)
    )

    # BUILD TMB object
    obj <- MakeADFun(data = data, parameters = params, hessian = TRUE,
      random = random, DLL = "SAMtool", silent = !verbose)

    #
    control <- list(iter.max = 5e3, eval.max = 1e4)
    opt_hess <- FALSE
    n_restart <- ifelse(opt_hess, 0, 1)
    
    # CALL minimizer TODO: CODE own optimize
    mod <- SAMtool:::optimize_TMB_model(obj, control, opt_hess, n_restart)

    # convergence flag
    conv <- ifelse(is.character(mod[[1]]), 0, 1)

    # EXTRACT outputs
    report <- obj$report(obj$env$last.par.best)
    
    out[[i]] <- list(F = report$F, Fstatus = report$F / report$FMSY,
      B = report$B[-ny], Bstatus = (report$B / report$BMSY)[-ny],
      Bdepletion = (report$B/report$K)[-ny],
      conv=conv, rps= unlist(report[c("FMSY", "BMSY", "K", "MSY")]))

    out[[i]]$data <- data
    out[[i]]$params <- params
    out[[i]]$mod <- mod
  }

  # AGGREGATE results
  ind <- FLQuants(lapply(setNames(nm=names(out[[1]][1:5])), function(x)
    empty %=% unlist(lapply(out, '[[', x))))

  rps <- Reduce(combine, lapply(out, function(x) FLPar(x$rps)))

  # WRITE tracking
  tracking['conv.est', ac(args$ay)] <- unlist(lapply(out, '[[', 'conv'))

  return(list(ind=ind, tracking=tracking, refpts=rps))
}
# }}}

# --- spict

# buildTMBinputSPICT {{{

# obsC, timeC, obsI, timeI

buildTMBinputSPICT <- function(catch, indices) {

  # PARSE and CHECK arguments
  ny <- dim(catch)[2]

  # C_hist

  dataC <- as.data.frame(catch)
  obsC <- dataC$data
  timeC <- dataC$year

  # obsI

  dataI <- lapply(indices, function(x) as.data.frame(x))
  obsI <- lapply(dataI, '[[', 'data')
  timeI <- lapply(dataI, '[[', 'year')

  # BUILD input list
  inp <- list(obsC=obsC, timeC=timeC, obsI=obsI, timeI=timeI)

  return(inp)
}
# }}}

# spict.sa {{{

spict.sa <- function(stk, idx, args, tracking) {
 
  # EXTRACT inputs
  catch <- catch(stk)
  indices <- lapply(idx, index)

  # TODO: DIMS
  ny <- dim(catch)[2]
  nsurvey <- length(indices)

  # TODO: CHECK inputs

  # OUTPUT
  empty <- catch %=% 0
  out <- vector(mode="list", length=args$it)

  # LOOP over iters
  for (i in seq(args$it)) {

    # EXTRACT single iter
    cab <- iter(catch, i)
    ind <- iter(indices, i)

    inp <- buildTMBinputSPICT(cab, ind)

    # SET simple options
    inp <- spict.simple(inp)

    # FIT spict
    fit <- fit.spict(inp)

    # convergence flag
    conv <- fit$opt$convergence

    # GET refpts
    rps <- unlist(c(fit$report[c("MSY", "Fmsy", "Bmsy")],
      K=fit$value[['K']]))

    # EXTRACT output indicators
    ind <- spict2ind(fit)

    out[[i]] <- list(ind=ind, conv=conv, rps=rps)
    
  }

  # COMBINE results

  inds <- lapply(out, "[[", "ind")
  ind <- FLQuants(lapply(setNames(nm=names(inds[[1]])), function(i)
    Reduce(combine, lapply(inds, "[[", i))))

  rps <- Reduce(combine, lapply(out, function(x) FLPar(x$rps)))

  # WRITE tracking
  tracking['conv.est', ac(args$ay)] <- unlist(lapply(out, '[[', 'conv'))

  return(list(ind=ind, tracking=tracking, refpts=rps))
}
# }}}

# spict.simple {{{

#' Function to simplify spict for a fast and robust as an MP
#'
#' - LIMITS obs CV, but can still esimate
#'
#' @param inp spict input time series
#' @param r.pr lognormal r prior, untransformed mean,log.sd, phase
#' @param bk.pr lognormal initial depletion prior
#' @param shape.pr lognormal prior for shape
#' @param pe lognormal prior for process error on biomass
#' @param oe lognormal prior for observation error on biomass
#' @param fdevs lognormal F penalty
#' @param ce lognorlam Catch penalty
#' @param dteuler time step resolution

spict.simple = function(inp, r.pr=c(0.3, 0.3, 1), bk.pr=c(0.9, 0.3, 1),
  shape.pr=c(2, 0.01, 1), oe=c(0.25, 0.1, 1), pe=c(0.07, 0.2, 1),
  fdevs=c(4, 0.5, 1), ce=c(0.1, 0.001, 1), dteuler=1) {

  # ASSIGN arguments to input
  inp$dteuler=dteuler

  # DEACTIVATE ratio of proc to obs
  inp$priors$logalpha <- c(0,0,0) 

  # DEACTIVATE catch to f devs
  inp$priors$logbeta <- c(0,0,0) 
  
  # SET bk prior
  inp$priors$logbkfrac <- 
    c(log(bk.pr[1])-bk.pr[2]^2/2,bk.pr[2],bk.pr[3])

  # SET obs error
  inp$priors$logsdi
  for(i in 1:length(inp$obsI)){
    inp$priors$logsdi[[i]] = c(log(oe[1])-0.5*oe[2]^2, oe[2], oe[3])  
  }

  # SET process error   
  inp$priors$logsdb <- c(log(pe[1])-0.5*pe[2]^2, pe[2], pe[3]) 

  # SET F devs: should be higher
  inp$priors$logsdf <- c(log(fdevs[1])-0.5*fdevs[2]^2, fdevs[[2]], fdevs[3])

  # SET catch error
  inp$priors$logsdc <- c(log(ce)[1]-0.5*ce[2]^2, ce[2], ce[3])

  # REDUCE shape CV
  inp$priors$logn <- c(log(shape.pr[1]),shape.pr[2],shape.pr[3])

  # SET r prior
  inp$priors$logr <- c(log(r.pr[1])-0.5*r.pr[2]^2,r.pr[2],r.pr[3])

  return(inp)
}

# }}}

# spict2FLQuant {{{

spict2FLQuant <- function(x,
  val=c("logB", "logFnotS", "logCpred", "logBBmsy", "logFFmsynotS")) {

  val <- match.arg(val)

  # logB
  if(val == "logB") {
    vec = x$par.random
  } else {
    vec = x$value
  }

  quant <- an(exp(vec[which(names(vec)==val)]))

  if(val == "logCpred") {
    season = x$inp$dtc[1]
  } else {
    season = round(1/x$inp$dteuler)
  }

  year <- c(rep(seq(x$inp$timerange[1], x$inp$timerange[2]), each=season))
  quant <- quant[1:length(year)]

  dat <- data.frame(age=1,
    year=c(rep(seq(x$inp$timerange[1], x$inp$timerange[2]),each=season)),
    season=c(rep(seq(season), x$inp$timerange[2] - x$inp$timerange[1]+1)),
    data=quant)

  if(val == "logCpred") {
    out <- seasonSums(as.FLQuant(dat))
  } else {
      out <- as.FLQuant(dat)[,,,season] 
      dimnames(out)$season <- "all"
  }
  return(out)
}
# }}}

# spict2ind {{{

spict2ind <- function(fit) {

  out <- list(
    B = spict2FLQuant(fit, val="logB"),
    F = spict2FLQuant(fit, val="logFnotS"),
    Bstatus = spict2FLQuant(fit, val="logBBmsy"),
    Fstatus = spict2FLQuant(fit, val="logFFmsynotS"),
    C = spict2FLQuant(fit, val="logCpred"))

  return(FLQuants(out))
}
# }}}

#    stk@refpts["B0"] = res$value["K"]/res$report$Bmsy
