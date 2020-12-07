### written by K. Garner, April 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

# ----------------------------------------------------------------------------------------------------
rm(list=ls())
# ----------------------------------------------------------------------------------------------------
### run analysis of sample size x effect size variability on the SRT data
# ----------------------------------------------------------------------------------------------------
# load packages and source function files
# -
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the location of this file
# uncomment the below and run if you need to install the packages
# install.packages("tidyverse")
# install.packages("wesanderson")
# install.packages("cowplot")
library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(ggridges)
library(lme4)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# load this guy if you have it already
load("SRT_sim_data.RData")

# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../total_of_313_subs_SRT_task_trial_level_data.csv", header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframes 
# ----------------------------------------------------------------------------------------------------

# Create a summary of the data for fixed fx modelling
min.RT <- 200 # in msec
sd.crit <- 2.5

ffx.dat <- dat %>% filter(Block.No > 2) %>%
              group_by(Subj.No, Block.No.Names) %>%
              filter(Accuracy == 1) %>%
              filter(RT.ms > min.RT) %>%
              filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
              summarise(RT=mean(RT.ms))

# now do the same for rfx modelling
# rfx.dat <- dat %>% filter(Block.No > 2) %>%
#               group_by(Subj.No, Block.No.Names) %>%
#               filter(Accuracy == 1) %>%
#               filter(RT.ms > min.RT) %>%
#               filter(RT.ms < (mean(RT.ms) + sd.crit*sd(RT.ms))) %>%
#               summarise(RT=mean(RT.ms))


# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1000# for each sample size, we will repeat our experiment n.perms times

# ----------------------------------------------------------------------------------------------------
# define variables for saving plots
# ----------------------------------------------------------------------------------------------------

plot.fname = "SRT.png"
rfx.plot.fname = "SRT_rfx.png"
width = 10 # in inches
height = 10

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from linear models, and cohen's d values, and save results to a list
# ----------------------------------------------------------------------------------------------------
subs = unique(dat$Subj.No)
sims = replicate(n.perms, lapply(sub.Ns, function(x) run.SRT.sim(data=ffx.dat, 
                                                                 dv="RT", 
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

# subs = unique(rfx.dat$Subj.No)
# rfx.sims = replicate(n.perms, lapply(sub.Ns, function(x) run.SRT.sim(data=rfx.dat, 
#                                                                      dv="RT", 
#                                                                      subs=subs,
#                                                                      N=x,
#                                                                      fx="rfx")), simplify = FALSE)
# 
# # ----------------------------------------------------------------------------------------------------
# # simplify and add to the sims.dat data.frame
# # ----------------------------------------------------------------------------------------------------
# 
# tmp = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, rfx.sims[[x]]))
# tmp = do.call(rbind, tmp)
# tmp$n <- as.factor(tmp$n)
# rfx = tmp %>% select(n, esub, eRes)
# tmp = tmp %>% select(-c('esub','eRes')) %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to = "value")
# tmp$measure <- as.factor(tmp$measure)
# tmp$model = "RFX"
# tmp$fx = "tt"
# sims.dat = rbind(sims.dat, tmp)
# rm(tmp)

# ----------------------------------------------------------------------------------------------------
# plot the outputs separately - then make 4 panels, top row = effect size, bottom row = p, left column = ffx, 
# right column = rfx
# ----------------------------------------------------------------------------------------------------

# first for d values
xlims = c(0,2)
sims.dat$fx = "tt"
ffx.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "FFX", ], xlims)

# now for p-values
xlims = c(0,.06)
ffx.p.p <- plt.ps(sims.dat[sims.dat$model=="FFX",], xlims)
# use cowplot to make a grid
p = plot_grid(ffx.d.p, ffx.p.p, labels=c('A', 'B'), label_size = 12, align="v")
# p # print out the plot so you can see it
p = p + ggsave(plot.fname, width = width/2, height = height, units="in")


# ----------------------------------------------------------------------------------------------------
# save the data of import to an RData file
# ----------------------------------------------------------------------------------------------------
save(sims.dat, file="SRT_sim_data.RData")
