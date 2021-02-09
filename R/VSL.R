### written by K. Garner, April 2020
### edited by Z. Nott, August 2020
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
library(ggridges)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# load this guy if you have it already
#load("VS_sim_data.RData")
set.seed(42) # testing diff seeds on output
# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../data/total_of_313_subs_VSL_task_trial_level_data.csv", header=TRUE)

# ----------------------------------------------------------------------------------------------------
# Create dataframe for t-test analysis
# ----------------------------------------------------------------------------------------------------

acc.dat <- dat %>% group_by(Subj.No) %>%
                   summarise(acc = mean(Accuracy))

# ----------------------------------------------------------------------------------------------------
# Create dataframe for permutations/prevalence analysis
# ----------------------------------------------------------------------------------------------------

# data frame contains TRUE ordering
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

sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1# for each sample size, we will repeat our experiment n.perms times
k = 1000 #for Monte Carlo simulations for prevalence stats (applies to both first level and second level perms)

# ----------------------------------------------------------------------------------------------------
# define variables for saving plots
# ----------------------------------------------------------------------------------------------------

plot.fname = "VSL.png"
width = 10 # in inches
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
flvl.perms <- run.mont.frst.lvl.over.subs(prev.dat, k) #P1 in 10.1016/j.neuroimage.2016.07.040 (label shuffle)
# now, over 100 experiments at each sample size, select the data, and then run the min stat procedure
prev.res <- replicate(n.perms, lapply(sub.Ns, function(x) run.prev.test(data=flvl.perms, 
                                                                            subs=subs,
                                                                            alpha=.05,
                                                                            N=x,
                                                                            k=k,
                                                                            Np=k)), simplify = FALSE)
prev.res <- do.call(rbind, do.call(rbind, prev.res))
# pivot longer and rename as p and d, and then add x-label to the plot below
prev.out <- prev.res %>% 
            pivot_longer(c('gamma', 'p'), names_to = "measure", values_to="value") %>%
            mutate(measure=recode(measure,
                                  'gamma' = 'd',
                                  'p' = 'p'))
prev.out$model = "Prevalence"
prev.out$n <- as.factor(prev.out$n)
# ----------------------------------------------------------------------------------------------------
# plot the outputs separately - then make 4 panels, top row = effect size, bottom row = p, left column = ffx, 
# right column = rfx
# ----------------------------------------------------------------------------------------------------

# first for d values
ylims = c(0, 2)
sims.dat$model <- "Null=.5"
sims.dat$fx <-  "t"
ffx.d.p <- plt.fx.sz(sims.dat, ylims)
ylims = c(0, 1)
prev.out$fx <- "min"
prev.d.p <- plt.fx.sz(prev.out, ylims) + xlab(expression(gamma))

# now for p-values
xlims=c(0,1)
ffx.p.p <- plt.ps(sims.dat, xlims)
prev.p.p <- plt.ps(prev.out, c(0,1))

# use cowplot to make a grid
p = plot_grid(ffx.d.p, prev.d.p, ffx.p.p, prev.p.p, labels=c('A', 'B', 'C', 'D'), label_size = 12, align="v")
# p # print out the plot so you can see it
p = p + ggsave(plot.fname, width = width, height = height, units="in")

# ----------------------------------------------------------------------------------------------------
# save the data of import to an RData file
# ----------------------------------------------------------------------------------------------------
save(sims.dat, prev.res, file="VS_sim_data.RData")