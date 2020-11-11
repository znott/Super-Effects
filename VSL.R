### written by K. Garner, April 2020
### edited by Z. Nott, August 2020
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
library(ggridges)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../total_of_313_subs_VSL_task_trial_level_data.csv", header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframe for t-test analysis
# ----------------------------------------------------------------------------------------------------

acc.dat <- dat %>% group_by(Subj.No) %>%
                   summarise(acc = mean(Accuracy))

# ----------------------------------------------------------------------------------------------------
# Create dataframe for permutations/prevalence analysis
# ----------------------------------------------------------------------------------------------------

prev.dat <- dat %>% select(Subj.No, Trial.No, Response, Target.Order, Accuracy) 
prev.dat$Response <- as.factor(prev.dat$Response)
prev.dat$Target.Order <- as.factor(prev.dat$Target.Order)

prev.dat <- prev.dat %>% mutate(Response = recode(Response,
                                      "122" = "Novel",
                                      "109" = "Repeat"),
                                Target.Order = recode(Target.Order,
                                      "1" = "Novel",
                                      "2" = "Repeat"))

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
# ----------------------------------------------------------------------------------------------------

sub.Ns = seq(23, 303, by = 10) 
n.perms =1000# for each sample size, we will repeat our experiment n.perms times
k = 1000 #for Monte Carlo simulations for prevalence stats

# ----------------------------------------------------------------------------------------------------
# define variables for saving plots
# ----------------------------------------------------------------------------------------------------

plot.fname = "VSL.png"
width = 5 # in inches
height = 10

# ----------------------------------------------------------------------------------------------------
# run simulations for t-test model, getting p values from t.tests, and cohen's d values, and save results to a list
# ----------------------------------------------------------------------------------------------------

subs = unique(acc.dat$Subj.No)
sims = replicate(n.perms, lapply(sub.Ns, function(x) run.os.t.test.sim(data=acc.dat, 
                                                                       dv="acc", 
                                                                       subs=subs,
                                                                       N=x)), simplify = FALSE)

# ----------------------------------------------------------------------------------------------------
# simplify the sims output to a dataframe and do a little wrangling to neaten it up
# ----------------------------------------------------------------------------------------------------

sims.dat = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, sims[[x]]))
sims.dat = do.call(rbind, sims.dat)
sims.dat$n <- as.factor(sims.dat$n)
sims.dat = sims.dat %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to="value")
sims.dat$measure <- as.factor(sims.dat$measure)

# ----------------------------------------------------------------------------------------------------
# minimum statistic approach
# ----------------------------------------------------------------------------------------------------

# first get the data from the first level perms
perm.dat <- run.mont.frst.lvl.over.subs(prev.dat, k)
prev.res <- replicate(n.perms, lapply(sub.Ns, function(x) run.scnd.lvl.mc(perm.dat, x, k)))
# pivot longer and rename as p and d, and then add x-label to the plot below
prev.out <- do.call(rbind, prev.res) %>% 
            pivot_longer(c('prev', 'p'), names_to = "measure", values_to="value") %>%
            mutate(measure=recode(measure,
                                  'prev' = 'd',
                                  'p' = 'p'))
prev.out$n <- as.factor(prev.out$N)
# ----------------------------------------------------------------------------------------------------
# plot the outputs separately - then make 4 panels, top row = effect size, bottom row = p, left column = ffx, 
# right column = rfx
# ----------------------------------------------------------------------------------------------------

# first for d values
ylims = c(0,3)
sims.dat$model <- "Null=.5"
ffx.d.p <- plt.fx.sz(sims.dat, ylims)
prev.d.p <- plt.fx.sz(prev.out, ylims) + xlab("prevalence")

# now for p-values
xlims=c(0,.06)
ffx.p.p <- plt.ps(sims.dat, xlims)
prev.p.p <- plt.ps(prev.out, c(0,1))

# use cowplot to make a grid
p = plot_grid(ffx.d.p, prev.d.p, ffx.p.p, prev.p.p, labels=c('A', 'B', 'C', 'D'), label_size = 12, align="v")
p # print out the plot so you can see it
p = p + ggsave(plot.fname, width = width, height = height, units="in")

# ----------------------------------------------------------------------------------------------------
# save the data of import to an RData file
# ----------------------------------------------------------------------------------------------------
save(sims.dat, rfx, file="VS_sim_data.RData")