### written by K. Garner, April 2020
### edited by Z. Nott, August 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

rm(list=ls())
### run analysis of sample size x effect size variability on the AB data
# ----------------------------------------------------------------------------------------------------
# load packages and source function files

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the location of this file
# uncomment the below and run if you need to install the packages
# install.packages("tidyverse")
# install.packages("wesanderson")
# install.packages("cowplot")
library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
library(car)
library(parallel)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting


# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
dat = read.csv("../data/total_of_313_subs_CC_task_trial_level_data.csv", header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframes 
# ----------------------------------------------------------------------------------------------------

# Create a summary of the data for fixed fx modelling
min.RT <- 200 # in msec
sd.crit <- 2.5

# really should just be 'dat', but I'm too committed to change it now
ffx.dat <- dat %>% mutate(Block.No = rep(c(1:12), each = 24, length(unique(dat$Subj.No)))) %>%
            group_by(Subj.No, Block.No, Trial.Type.Name) %>%
            filter(Accuracy == 1) %>%
            filter(RT.ms > min.RT) %>%
            filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
            summarise(RT=mean(RT.ms))

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1000# for each sample size, we will repeat our experiment n.perms times
cores = 20

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from t.tests, and cohen's d values, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs  <- unique(ffx.dat$Subj.No)
lapply(sub.Ns, function(x) run.outer(in.data=ffx.dat, subs=subs, N=x, k=n.perms, j=n.perms, cores=cores, ffx.f=get.ps.CC, rfx.f=run.lme.4.cc, fstem="CC_N-%d_parent-%d.RData"))

# ----------------------------------------------------------------------------------------------------
# attain densities for each subject N, across all outer samples
# ----------------------------------------------------------------------------------------------------
dens.across.N(fstem="CC_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=-50, max=0, spacer=1000, dv="p", savekey="CC")
dens.across.N(fstem="CC_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=0.5, spacer=1000, dv="d", savekey="CC")
dens.across.N(fstem="CC_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=400, spacer=1000, dv="esub", savekey="CC")
dens.across.N(fstem="CC_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=400, spacer=1000, dv="eRes", savekey="CC")

quit()

