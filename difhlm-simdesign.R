#06mar2018

rm(list=ls())

sessionInfo()

library(SimDesign)
library(dplyr)
library(CTT)
library(lme4)
source("difhlm-simdesign-functions.R")

# set up item parameters and simulation conditions
ip   <- makeItemParms(startparms = "item-parms-mcas-2015-math-gr8.csv")
pars <- list(iparms=ip, fititem=c(25)) # fititem indicates which items to test for DIF on each replication
reps <- 200

Design <- expand.grid(
  icc = c(0.05, 0.3),
  J = c(25, 100),
  N = c(10, 40),
  impact = c(0, 1),
  dif = c("no", "with"),
  difvar = c("no", "with"),
  bval = c("low", "med", "high")
)

set.seed(20171214)
use_seeds <- ifelse(runif(nrow(Design))>0.5, -1, 1)*round(runif(nrow(Design))*1e8)
use_seeds

stopifnot(length(unique(use_seeds)) == nrow(Design))

# Run
result <- runSimulation(
  Design,
  replications  = reps,
  parallel      = TRUE,
  fixed_objects = pars,
  generate      = Generate,
  analyse       = Analyze,
  summarise     = Summarize,
  as.factor     = FALSE,
  save_results  = TRUE,
  save_details  = list(
    save_results_dirname = "difhlm-simdesign-results",
    save_seeds_dirname   = "difhlm-simdesign-seeds"),
  save       = TRUE,
  filename   = "difhlm-simdesign-results",
  save_seeds = TRUE,
  seed = use_seeds)

sessionInfo()

