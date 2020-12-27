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
library(lme4)
library(ggridges)
library(car)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

#$ load a previous state if you have it
load("CC_sim_data.RData")

# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
dat = read.csv("../total_of_313_subs_CC_task_trial_level_data.csv", header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframes 
# ----------------------------------------------------------------------------------------------------

# Create a summary of the data for fixed fx modelling
min.RT <- 200 # in msec
sd.crit <- 2.5

ffx.dat <- dat %>% mutate(Block.No = rep(c(1:12), each = 24, length(unique(dat$Subj.No)))) %>%
            group_by(Subj.No, Block.No, Trial.Type.Name) %>%
            filter(Accuracy == 1) %>%
            filter(RT.ms > min.RT) %>%
            filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
            summarise(RT=mean(RT.ms))

# now do the same for rfx modelling
rfx.dat <- dat %>% mutate(Block.No = rep(c(1:12), each = 24, length(unique(dat$Subj.No)))) %>%
            group_by(Subj.No, Task.Order, Experimenter, Block.No, Trial.Type.Name) %>%
            filter(Accuracy == 1) %>%
            filter(RT.ms > min.RT) %>%
            filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
            summarise(RT=mean(RT.ms)) 

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1000# for each sample size, we will repeat our experiment n.perms times

# ----------------------------------------------------------------------------------------------------
# define variables for saving plots
# ----------------------------------------------------------------------------------------------------

plot.fname = "CC"
rfx.plot.fname = "CC_rfx.png"
width = 10 # in inches
height = 10

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from t.tests, and cohen's d values, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs = unique(ffx.dat$Subj.No)
sims = replicate(n.perms, lapply(sub.Ns, function(x) run.aov.CC.sim(data=ffx.dat, subs=subs,
                                                                    N=x, dv="RT",fx="ffx")), simplify = FALSE)

# ----------------------------------------------------------------------------------------------------
# simplify the sims output to a dataframe and do a little wrangling to neaten it up
# ----------------------------------------------------------------------------------------------------
sims.dat = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, sims[[x]]))
sims.dat = do.call(rbind, sims.dat)
sims.dat = sims.dat %>% select(-c('esub', 'eRes'))
sims.dat$n <- as.factor(sims.dat$n)
sims.dat = sims.dat %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to="value")
sims.dat$measure <- as.factor(sims.dat$measure)
sims.dat$model <- "FFX"

# ----------------------------------------------------------------------------------------------------
# run simulations for rfx **interaction** models, getting p values and d for rfx, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs = unique(rfx.dat$Subj.No)
rfx.sims = replicate(n.perms, lapply(sub.Ns, function(x) run.aov.CC.sim(data=rfx.dat, 
                                                                        dv="RT", 
                                                                        subs=subs,
                                                                        N=x,
                                                                        fx="rfx",
                                                                        efx="int")), simplify = FALSE)

# ----------------------------------------------------------------------------------------------------
# simplify and add to the sims.dat data.frame
# ----------------------------------------------------------------------------------------------------

tmp = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, rfx.sims[[x]]))
tmp = do.call(rbind, tmp)
tmp$n <- as.factor(tmp$n)
rfx = tmp %>% select(n, esub, eRes)
rfx$model = "int"
tmp = tmp %>% select(-c('esub','eRes')) %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to = "value")
tmp$measure <- as.factor(tmp$measure)
tmp$model = "RFX"
sims.dat = rbind(sims.dat, tmp)
rm(tmp)

# ----------------------------------------------------------------------------------------------------
# run simulations for rfx **ME** models, getting p values and d squares for rfx, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs = unique(rfx.dat$Subj.No)
rfx.sims = replicate(n.perms, lapply(sub.Ns, function(x) run.aov.CC.sim(data=rfx.dat, 
                                                                        dv="RT", 
                                                                        subs=subs,
                                                                        N=x,
                                                                        fx="rfx",
                                                                        efx="me")), simplify = FALSE)

# ----------------------------------------------------------------------------------------------------
# simplify and add to the sims.dat data.frame
# ----------------------------------------------------------------------------------------------------

tmp = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, rfx.sims[[x]]))
tmp = do.call(rbind, tmp)
tmp$n <- as.factor(tmp$n)
tmprfx = tmp %>% select(n, esub, eRes)
tmprfx$model = "me"
rfx = rbind(rfx, tmprfx)
tmp = tmp %>% select(-c('esub','eRes')) %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to = "value")
tmp$measure <- as.factor(tmp$measure)
tmp$model = "RFX"
sims.dat = rbind(sims.dat, tmp)
rm(tmp)
rm(tmprfx)


# ----------------------------------------------------------------------------------------------------
# plot the outputs separately - then make 4 panels, top row = effect size, bottom row = p, left column = ffx, 
# right column = rfx
# ----------------------------------------------------------------------------------------------------

# first for d values
ylims = c(0,0.75)
ffxme.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "FFX" & sims.dat$fx == "me", ], ylims)
ffxint.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "FFX" & sims.dat$fx == "int", ], ylims)
rfxint.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "RFX" & sims.dat$fx == "me", ], ylims)
rfxme.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "RFX" & sims.dat$fx == "int", ], ylims)

# now for p-values
ffxme.p.p <- plt.ps(sims.dat[sims.dat$model=="FFX" & sims.dat$fx == "me",], c(0,1))
ffxint.p.p <- plt.ps(sims.dat[sims.dat$model=="FFX" & sims.dat$fx == "int",], c(0,1))
rfxint.p.p <- plt.ps(sims.dat[sims.dat$model=="RFX" & sims.dat$fx == "me",], c(0,1))
rfxme.p.p <- plt.ps(sims.dat[sims.dat$model=="RFX"& sims.dat$fx == "int",], c(0,1))

# use cowplot to make a grid
dp = plot_grid(ffxme.d.p, ffxint.d.p, rfxint.d.p, rfxme.d.p, labels=c('A', 'B', 'C', 'D'), label_size = 12, align="v")
# # print out the plot so you can see it
dp = dp + ggsave(paste(plot.fname, "d.png", sep="_"), width = width, height = height, units="in")

pp = plot_grid(ffxme.p.p, ffxint.p.p, rfxint.p.p, rfxme.p.p, labels=c('A', 'B', 'C', 'D'), label_size = 12, align="v")
pp = pp + ggsave(paste(plot.fname, "p.png", sep="_"), width = width, height = height, units="in")

# ----------------------------------------------------------------------------------------------------
# now a raincloud plot of the sources of randomness in the model
# ----------------------------------------------------------------------------------------------------

rfx.p <- plt.rfx(rfx, c(0, 50000))
rfx.p = rfx.p + ggsave(rfx.plot.fname, width = width, height = height, units="in")

# ----------------------------------------------------------------------------------------------------
# save the data of import to an RData file
# ----------------------------------------------------------------------------------------------------
save(sims.dat, rfx, file="CC_sim_data.RData")


