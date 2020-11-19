### written by K. Garner, April 2020
### edited by Z. Nott, July 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

rm(list=ls())
### run analysis of sample size x effect size variability on the SRT data
# ----------------------------------------------------------------------------------------------------
# load packages and source function files

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the location of this file
# uncomment the below and run if you need to install the packages
# install.packages("tidyverse")
# install.packages("wesanderson")
# install.packages("cowplot")
library(tidyverse) # for data wrangling
library(wesanderson) # palette for some sweet figure colours
library(cowplot)
library(lme4)
library(ggridges)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../total_of_313_subs_SingDual_task_trial_level_data.csv", header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframes 
# ----------------------------------------------------------------------------------------------------
# Create a summary of the data for ffx and rfx modelling
min.RT <- .200 # in sec
sd.crit <- 2.5

rfx.dat <- dat %>% filter(Overall.Accuracy == 1) %>%
                   select(c('Subj.No', 'Task.Order', 'Experimenter', 'Trial.Type.Name', 'Task.1.RT.Sec', 'Task.2.RT.Sec'))  %>%
                   pivot_longer(c('Task.1.RT.Sec', 'Task.2.RT.Sec'), names_to = "task", values_to="RT") %>%
                   drop_na() 
rfx.dat$task[rfx.dat$Trial.Type.Name == 'single_auditory'] = 'sound'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'single_visual'] = 'vis'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'dual_task' & rfx.dat$task == 'Task.1.RT.Sec'] = 'vis'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'dual_task' & rfx.dat$task == 'Task.2.RT.Sec'] = 'sound'

rfx.dat <- rfx.dat %>% mutate(trialtype = fct_recode(Trial.Type.Name,
                                                                    'single' = 'single_auditory',
                                                                    'single' = 'single_visual',
                                                                    'dual' = 'dual_task')) %>%
                       select(-Trial.Type.Name) %>%
                       group_by(Subj.No, Task.Order, Experimenter, task, trialtype) %>%
                       filter(RT > min.RT) %>%
                       filter(RT < (mean(RT)+sd.crit*sd(RT))) %>%
                       summarise(RT = mean(RT))

ffx.dat <- rfx.dat %>% ungroup(Experimenter, Task.Order) %>% select(-c('Task.Order', 'Experimenter'))

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1000# for each sample size, we will repeat our experiment n.perms times

# ----------------------------------------------------------------------------------------------------
# define variables for saving plots
# ----------------------------------------------------------------------------------------------------

plot.fname = "SD.png"
rfx.plot.fname = "SD_rfx.png"
width = 10 # in inches
height = 10

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from linear models, and cohen's d values, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs = unique(ffx.dat$Subj.No)
sims = replicate(n.perms, lapply(sub.Ns, function(x) run.SD.sim(data=ffx.dat,
                                                                     dv="RT",
                                                                     subs=subs,
                                                                     N=x,
                                                                     fx='ffx')), simplify = FALSE)

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
rfx.sims = replicate(n.perms, lapply(sub.Ns, function(x) run.SD.sim(data=rfx.dat, 
                                                                     dv="RT", 
                                                                     subs=subs,
                                                                     N=x,
                                                                     fx="rfx")), simplify = FALSE)

# ----------------------------------------------------------------------------------------------------
# simplify and add to the sims.dat data.frame, make a separate dataframe of the random fx info
# ----------------------------------------------------------------------------------------------------

tmp = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, rfx.sims[[x]]))
tmp = do.call(rbind, tmp)
tmp$n <- as.factor(tmp$n)
rfx = tmp %>% select(n, esub, eRes)
tmp = tmp %>% select(-c('esub','eRes')) %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to = "value")
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
ffx.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "FFX", ], ylims)
rfx.d.p <- plt.fx.sz(sims.dat[sims.dat$model == "RFX", ], ylims)

# now for p-values
ffx.p.p <- plt.ps(sims.dat[sims.dat$model=="FFX",])
rfx.p.p <- plt.ps(sims.dat[sims.dat$model=="RFX",])

# use cowplot to make a grid
p = plot_grid(ffx.d.p, rfx.d.p, ffx.p.p, rfx.p.p, labels=c('A', 'B', 'C', 'D'), label_size = 12, align="v")
# p # print out the plot so you can see it
p = p + ggsave(plot.fname, width = width, height = height, units="in")


# ----------------------------------------------------------------------------------------------------
# now a raincloud plot of the sources of randomness in the model
# ----------------------------------------------------------------------------------------------------

rfx$esub = rfx$esub*1000 # putting in ms
rfx$eRes = rfx$eRes*1000
rfx.p <- plot.rfx(rfx, c(0, 25))
rfx.p = rfx.p + ggsave(rfx.plot.fname, width = width/2, height = height/2, units="in")


# ----------------------------------------------------------------------------------------------------
# save the data of import to an RData file
# ----------------------------------------------------------------------------------------------------
save(sims.dat, rfx, file="SD_sim_data.RData")

