### written by K. Garner, April 2020
### edited by Z. Nott, July 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

# ----------------------------------------------------------------------------------------------------
rm(list=ls())
# ----------------------------------------------------------------------------------------------------
### run analysis of sample size x effect size variability on the AB data
# ----------------------------------------------------------------------------------------------------
# load packages and source function files
# ----------------------------------------------------------------------------------------------------

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the location of this file
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

#$ load a previous state if you have it
# load("AB_sim_data.RData")

# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../data/total_of_313_subs_AB_task_trial_level_data.csv", header=TRUE)
dat$Task.Order <- as.factor(dat$Task.Order)
dat$Experimenter <- as.factor(dat$Experimenter)
dat$T1.Stimulus.Type <- as.factor(dat$T1.Stimulus.Type)
dat$T2.Stimulus.Type <- as.factor(dat$T2.Stimulus.Type)

# ----------------------------------------------------------------------------------------------------
# Create dataframes 
# ----------------------------------------------------------------------------------------------------

# Create a summary of the data for fixed fx modelling
ffx.dat <- dat %>% group_by(Subj.No, Trial.Type.Name) %>%
                   summarise(T1=mean(T1.Accuracy),
                             T2gT1=mean(T2T1.Accuracy))

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------
sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1000# for each sample size, we will repeat our experiment n.perms^2 times 
cores = 4

# ----------------------------------------------------------------------------------------------------
# define variables for saving plots
# ----------------------------------------------------------------------------------------------------

plot.fname = "AB.png"
rfx.plot.fname = "AB_rfx.png"
width = 10 # in inches
height = 10

# ----------------------------------------------------------------------------------------------------
# run simulations for ffx & rfx models, getting p values and partial eta squares, and save results to a list
# ----------------------------------------------------------------------------------------------------
subs  <- unique(ffx.dat$Subj.No)
# implement with fstem, and functions
lapply(sub.Ns, function(x) run.outer(in.data=ffx.dat, subs=subs, N=x, k=n.perms, j=n.perms, cores=cores, ffx.f=get.ps.aov.AB, rfx.f=run.lme.4.AB, fstem="AB_N-%d_parent-%d.RData"))

# ----------------------------------------------------------------------------------------------------
# attain densities for each subject N, across all outer samples
# ----------------------------------------------------------------------------------------------------
dens.across.N(fstem="AB_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=-800, max=0, spacer=1000, dv="p", savekey="AB")
dens.across.N(fstem="AB_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=3, spacer=1000, dv="d", savekey="AB")
dens.across.N(fstem="AB_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=1, spacer=1000, dv="esub", savekey="AB")
dens.across.N(fstem="AB_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=1, spacer=1000, dv="eRes", savekey="AB")


# ----------------------------------------------------------------------------------------------------
# plot the outputs separately - then make 4 panels, top row = effect size, bottom row = p, left column = ffx, 
# right column = rfx
# ----------------------------------------------------------------------------------------------------

# first for d values
# ylims = c(0,3)
# ffx.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "FFX", ], c(0,2))
# names(sims.dat)[names(sims.dat)=="rfx"] = "fx"
# rfx.sub.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "RFX" & sims.dat$fx == "sub", ], c(1,3))
# rfx.stim.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "RFX" & sims.dat$fx == "stim", ], c(0,2))
# # now for p-values
# xlims=c(0,0.4)
# ffx.p.p <- plt.ps(sims.dat[sims.dat$model=="FFX",], xlims)
# rfxsub.p.p <- plt.ps(sims.dat[sims.dat$model=="RFX" & sims.dat$fx == "sub",], c(0, .4)) + geom_density_ridges()
# rfxstim.p.p <- plt.ps(sims.dat[sims.dat$model=="RFX" & sims.dat$fx == "stim",], c(0, .4))
# 
# # use cowplot to make a grid
# 
# p = plot_grid(ffx.d.p, rfx.sub.d.p, rfx.stim.d.p, ffx.p.p, rfxsub.p.p, rfxstim.p.p, labels=c('A', 'B', 'C', 'D', 'E', 'F'), label_size = 12, align="v")
# #p # print out the plot so you can see it
# p = p + ggsave(plot.fname, width = width, height = height, units="in")
# 
# # ----------------------------------------------------------------------------------------------------
# # now a raincloud plot of the sources of randomness in the model
# # ----------------------------------------------------------------------------------------------------
# 
# rfx$model <- NULL
# names(rfx)[names(rfx) == "rfx"] = "model"
# rfx.p <- plt.rfx(rfx, c(0, 0.2)) + facet_wrap(~model*rfx, ncol=2)
# rfx.p = rfx.p + ggsave(rfx.plot.fname, width = width/2, height = height/2, units="in")
# 
# # ----------------------------------------------------------------------------------------------------
# # save the data of import to an RData file
# # ----------------------------------------------------------------------------------------------------
# save(sims.dat, rfx, file="AB_sim_data.RData")
