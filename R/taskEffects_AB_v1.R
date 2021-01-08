### written by K. Garner, April 2020
### edited by Z. Nott, July 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

rm(list=ls())
### run analysis of sample size x effect size variability on the AB data
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
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting

# ----------------------------------------------------------------------------------------------------
# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
dat = read.csv("EFIL_313_subs_AB_edit1.csv", header=TRUE)
dat$Subjects = as.factor(dat$Subjects)
dat = dat %>% pivot_longer(c('shortLag', 'longLag'), names_to = "block", values_to="RT") # move to a tidy frame
dat$block = as.factor(dat$block)
dat = dat %>% mutate(block = fct_recode(block, # rename the block factor and replace into dataframe
                                        "short" = "shortLag", 
                                        "long" = "longLag"))

# ----------------------------------------------------------------------------------------------------
# define levels for simulations
sub.Ns = seq(23, 303, by = 10) 
n.perms =1000# for each sample size, we will repeat our experiment n.perms times

# ----------------------------------------------------------------------------------------------------
# define variables for saving plots
plot.fname = "AB.png"
width = 5 # in inches
height = 10

# ----------------------------------------------------------------------------------------------------
# run simulations, getting p values from t.tests, and cohen's d values, and save results to a list
subs = unique(dat$Subjects)
sims = replicate(n.perms, lapply(sub.Ns, function(x) run.t.test.sim(data=dat, iv="block", 
                                                                    dv="RT", x="short", 
                                                                    y="long", subs=subs,
                                                                    N=x)), simplify = FALSE)

# ----------------------------------------------------------------------------------------------------
# simplify the sims output to a dataframe and do a little wrangling to neaten it up
sims.dat = lapply(seq(1,n.perms,by=1), function(x) do.call(rbind, sims[[x]]))
sims.dat = do.call(rbind, sims.dat)
sims.dat$n <- as.factor(sims.dat$n)
sims.dat = sims.dat %>% pivot_longer(c('p', 'd'), names_to = "measure", values_to="value")
sims.dat$measure <- as.factor(sims.dat$measure)

# ----------------------------------------------------------------------------------------------------
# plot the outputs separately - then make 2 panels, 1 with sample size x p-value, 1 with sample size x effect size

# first for d values
d.p <- ggplot(sims.dat[sims.dat$measure == "d", ], aes(x=n, y=value, fill = n, colour = n)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim =
                     TRUE) +
  geom_boxplot(aes(x = as.numeric(n)+0.25, y = value), outlier.shape = NA,
               alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('d') + xlab('N') + theme_cowplot() + 
  guides(fill = FALSE, colour = FALSE) +
  coord_flip() +           
  theme(axis.title.x = element_text(face = "italic"))

# now for p-values
p.p <- ggplot(sims.dat[sims.dat$measure == "p", ], aes(x=n, y=value, fill = n, colour = n)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim =
                     TRUE) +
  geom_boxplot(aes(x = as.numeric(n)+0.25, y = value), outlier.shape = NA,
               alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('p') + xlab('') + theme_cowplot() + ylim(c(0,1)) +
  geom_hline(aes(yintercept=.05), linetype="dashed") +
  guides(fill = FALSE, colour = FALSE) +
  coord_flip() +           
  theme(axis.title.x = element_text(face = "italic"),
        axis.text.y = element_blank())


# use cowplot to make a grid
p = plot_grid(d.p, p.p, labels=c('A', 'B'), label_size = 12)
p # print out the plot so you can see it
p = p + ggsave(plot.fname, width = width, height = height, units="in")




