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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the location of this file
# uncomment the below and run if you need to install the packages
# install.packages("tidyverse")
# install.packages("wesanderson")
# install.packages("cowplot")
library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4) # for mixed effects modelling
library(ggridges)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

#$ load a previous state if you have it
load("AB_sim_data.RData")

# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../total_of_313_subs_AB_task_trial_level_data.csv", header=TRUE)
dat$Task.Order <- as.factor(dat$Task.Order)
dat$Experimenter <- as.factor(dat$Experimenter)

# ----------------------------------------------------------------------------------------------------
# Summarise random effects for intel for lme modelling 
# ----------------------------------------------------------------------------------------------------

rfx.intel.task.order <- dat %>% group_by(Task.Order) %>%
                             summarise(count=length(unique(Subj.No)))
rfx.intel.experimenter <- dat %>% group_by(Experimenter) %>%
                             summarise(count=length(unique(Subj.No)))

# ----------------------------------------------------------------------------------------------------
# Create dataframes 
# ----------------------------------------------------------------------------------------------------

# Create a summary of the data for fixed fx modelling
ffx.dat <- dat %>% group_by(Subj.No, Trial.Type.Name) %>%
                   summarise(T1=mean(T1.Accuracy),
                             T2gT1=mean(T2T1.Accuracy))

# now do the same for rfx modelling
rfx.dat <- dat %>% group_by(Subj.No, Task.Order, Experimenter, Trial.Type.Name) %>%
                   summarise(T1=mean(T1.Accuracy),
                   T2gT1=mean(T2T1.Accuracy))

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1000# for each sample size, we will repeat our experiment n.perms times

# ----------------------------------------------------------------------------------------------------
# define variables for saving plots
# ----------------------------------------------------------------------------------------------------

plot.fname = "AB.png"
rfx.plot.fname = "AB_rfx.png"
width = 10 # in inches
height = 10

# ----------------------------------------------------------------------------------------------------
# run simulations for ffx models, getting p values and partial eta squares, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs = unique(ffx.dat$Subj.No)
sims = replicate(n.perms, lapply(sub.Ns, function(x) run.aov.AB.sim(data=ffx.dat, 
                                                                    dv="T2gT1", 
                                                                    subs=subs,
                                                                    N=x,
                                                                    fx="ffx")), simplify = FALSE)

# ----------------------------------------------------------------------------------------------------
# simplify the sims output to a dataframe and do a little wrangling to neaten it up
# ----------------------------------------------------------------------------------------------------

sims.dat = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, sims[[x]]))
sims.dat = do.call(rbind, sims.dat)
sims.dat = sims.dat %>% select(-c('esub', 'eRes'))
sims.dat$n <- as.factor(sims.dat$n)
sims.dat = sims.dat %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to="value")
sims.dat$measure <- as.factor(sims.dat$measure)
sims.dat$model = "FFX"

# ----------------------------------------------------------------------------------------------------
# run simulations for rfx models, getting p values and partial eta squares for ffx, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs = unique(rfx.dat$Subj.No)
rfx.sims = replicate(n.perms, lapply(sub.Ns, function(x) run.aov.AB.sim(data=rfx.dat, 
                                                                    dv="T2gT1", 
                                                                    subs=subs,
                                                                    N=x,
                                                                    fx="rfx")), simplify = FALSE)

# ----------------------------------------------------------------------------------------------------
# simplify and add to the sims.dat data.frame
# ----------------------------------------------------------------------------------------------------

tmp = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, rfx.sims[[x]]))
tmp = do.call(rbind, tmp)
tmp$n <- as.factor(tmp$n)
rfx = tmp %>% select(n, esub, eRes)
tmp = tmp %>% select(-c('esub', 'eRes')) %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to = "value")
tmp$measure <- as.factor(tmp$measure)
tmp$model = "RFX"
sims.dat = rbind(sims.dat, tmp)
rm(tmp)

# ----------------------------------------------------------------------------------------------------
# plot the outputs separately - then make 4 panels, top row = effect size, bottom row = p, left column = ffx, 
# right column = rfx
# ----------------------------------------------------------------------------------------------------

# first for d values
ylims = c(0,3)
ffx.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "FFX", ], c(0,2))
rfx.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "RFX", ], c(1,3))

# now for p-values
xlims=c(0,1)
ffx.p.p <- plt.ps(sims.dat[sims.dat$model=="FFX",], xlims)
rfx.p.p <- plt.ps(sims.dat[sims.dat$model=="RFX",], xlims)

# use cowplot to make a grid
p = plot_grid(ffx.d.p, rfx.d.p, ffx.p.p, rfx.p.p, labels=c('A', 'B', 'C', 'D'), label_size = 12, align="v")
#p # print out the plot so you can see it
p = p + ggsave(plot.fname, width = width, height = height, units="in")

# ----------------------------------------------------------------------------------------------------
# now a raincloud plot of the sources of randomness in the model
# ----------------------------------------------------------------------------------------------------

rfx.p <- plt.rfx(rfx, c(0, .025))
rfx.p = rfx.p + ggsave(rfx.plot.fname, width = width/2, height = height/2, units="in")

# ----------------------------------------------------------------------------------------------------
# save the data of import to an RData file
# ----------------------------------------------------------------------------------------------------
save(sims.dat, rfx, file="AB_sim_data.RData")
