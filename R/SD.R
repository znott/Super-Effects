### written by K. Garner, April 2020
### edited by Z. Nott, July 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

rm(list=ls())
### run analysis of sample size x effect size variability on the SRT data
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
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../data/total_of_313_subs_SingDual_task_trial_level_data.csv", header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframes 
# ----------------------------------------------------------------------------------------------------
# Create a summary of the data for ffx and rfx modelling
min.RT <- .200 # in sec
sd.crit <- 2.5

rfx.dat <- dat %>% filter(Overall.Accuracy == 1) %>%
                   select(c('Subj.No', 'Trial.Type.Name', 'Task.1.RT.Sec', 'Task.2.RT.Sec', 'Task.1.Response', 'Task.2.Response'))  %>%
                   pivot_longer(c('Task.1.RT.Sec', 'Task.2.RT.Sec'), names_to = "task", values_to="RT") %>%
                   drop_na()
rfx.dat$task[rfx.dat$Trial.Type.Name == 'single_auditory'] = 'sound'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'single_visual'] = 'vis'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'dual_task' & rfx.dat$task == 'Task.1.RT.Sec'] = 'vis'
rfx.dat$task[rfx.dat$Trial.Type.Name == 'dual_task' & rfx.dat$task == 'Task.2.RT.Sec'] = 'sound'

rfx.dat <- rfx.dat %>% mutate(trialtype = fct_recode(Trial.Type.Name,
                                                                    'single' = 'single_auditory',
                                                                    'single' = 'single_visual',
                                                                    'dual' = 'dual_task')) 
rfx.dat$task.stim <- NA
rfx.dat$task.stim[rfx.dat$trialtype == "single"] = rfx.dat$Task.1.Response[rfx.dat$trialtype == "single"]
rfx.dat$task.stim[rfx.dat$trialtype == "dual" & rfx.dat$task == "vis"] = rfx.dat$Task.1.Response[rfx.dat$trialtype == "dual" & rfx.dat$task == "vis"]
rfx.dat$task.stim[rfx.dat$trialtype == "dual" & rfx.dat$task == "sound"] = rfx.dat$Task.2.Response[rfx.dat$trialtype == "dual" & rfx.dat$task == "sound"]

rfx.dat <- rfx.dat %>% select(-c("Trial.Type.Name", "Task.1.Response", "Task.2.Response")) %>%
                       group_by(Subj.No, task, trialtype, task.stim) %>%
                       filter(RT > min.RT) %>%
                       filter(RT < (mean(RT)+sd.crit*sd(RT))) %>%
                       summarise(RT = mean(RT))

ffx.dat <- rfx.dat %>% ungroup(task.stim) %>% select(-c('task.stim'))

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1000# for each sample size, we will repeat our experiment n.perms times
cores = 20

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from linear models, and cohen's d values, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs  <- unique(ffx.dat$Subj.No)
lapply(sub.Ns, function(x) run.outer(in.data=ffx.dat, subs=subs, N=x, k=n.perms, j=n.perms, cores=cores, ffx.f=get.ps.SD, rfx.f=run.lme.4.SD, fstem="SD_N-%d_parent-%d.RData"))

# ----------------------------------------------------------------------------------------------------
# attain densities for each subject N, across all outer samples
# ----------------------------------------------------------------------------------------------------
dens.across.N(fstem="SD_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=-800, max=0, spacer=1000, dv="p", savekey="SD")
dens.across.N(fstem="SD_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=3, spacer=1000, dv="d", savekey="SD")
dens.across.N(fstem="SD_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=400, spacer=1000, dv="esub", savekey="SD")
dens.across.N(fstem="SD_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=400, spacer=1000, dv="eRes", savekey="SD")

# ----------------------------------------------------------------------------------------------------
# get outta here
# ----------------------------------------------------------------------------------------------------
quit()
