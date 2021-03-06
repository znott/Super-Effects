### written by K. Garner, April 2020
### edited by Z. Nott, Oct 2020
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
library(lme4)
source("efilids_functions.R") # custom functions written for this project
source("R_rainclouds.R") # functions for plotting


# ----------------------------------------------------------------------------------------------------
# load data and wrangle 
# ----------------------------------------------------------------------------------------------------
dat = read.csv("../total_of_313_subs_CC_task_trial_level_data.csv", header=TRUE)


##### RT_ms: Trial Type x Task Order
##### --------------------------------------------------------
task.order.cc <- dat %>% group_by(Subj.No, Task.Order, Block.No, Trial.Type) %>% 
                                  filter(RT.ms>.2) %>%
                                  filter(RT.ms < (mean(RT.ms) + 2.5*sd(RT.ms))) %>%
                                  summarise(mean=mean(RT.ms))
task.order.cc$Trial.Type <- factor(task.order.cc$Trial.Type)
task.order.cc$Task.Order <- factor(task.order.cc$Task.Order)
#task.order.cc$Subj.No <- factor(task.order.cc$Subj.No)
task.order.cc$Block.No <- factor(task.order.cc$Block.No)

t.o <- task.order.cc %>% ggplot(aes(x=Block.No, y=mean, group=Trial.Type)) +
  geom_line() +  ylim(500,2000) +
  facet_wrap(~Task.Order)
t.o

# define variables for saving plots
plot.fname = "CC_FX_Task.Order.png"
width = 30 # in inches
height = 10

#save output
ggsave(plot.fname, width = width, height = height, units="in", limitsize = FALSE)


##### RT_ms: Trial Type x Experimenter
##### --------------------------------------------------------
experimenter.cc <- dat %>% group_by(Experimenter, Block.No, Trial.Type) %>% 
                           filter(RT.ms>.2) %>%
                           filter(RT.ms < (mean(RT.ms) + 2.5*sd(RT.ms))) %>%
                           summarise(mean=mean(RT.ms))

experimenter.cc$Trial.Type <- factor(experimenter.cc$Trial.Type)
#experimenter.cc$Subj.No <- factor(experimenter.cc$Subj.No)
experimenter.cc$Block.No <- factor(experimenter.cc$Block.No)

ex <- experimenter.cc %>% ggplot(aes(x=Block.No, y=mean, group=Trial.Type, fill=Trial.Type, col=Trial.Type)) +
      geom_line() + ylim(500,1500) +
      facet_wrap(~Experimenter)
ex

# define variables for saving plots
plot.fname = "CC_FX_Experimenter.png"
width = 30 # in inches
height = 10

#save output
ggsave(plot.fname, width = width, height = height, units="in", limitsize = FALSE)


### RT_ms: Trial Type x Stimulus Type
##### --------------------------------------------------------
stimulus.cc <- dat %>% group_by(Subj.No, Stimulus.Type, Block.No, Trial.Type, RT.ms) %>% summarise(mean=mean(RT.ms))
stimulus.cc$Trial.Type <- factor(stimulus.cc$Trial.Type)
stimulus.cc$Stimulus.Type <- factor(stimulus.cc$Stimulus.Type)
stimulus.cc$Subj.No <- factor(stimulus.cc$Subj.No)
stimulus.cc$Block.No <- factor(stimulus.cc$Block.No)

st <- stimulus.cc %>% ggplot(aes(x=Stimulus.Type, y=RT.ms, fill=Trial.Type)) + ylim(0,3000) +
  geom_boxplot() +
  facet_wrap(~Block.No)
st

# define variables for saving plots
plot.fname = "CC_FX_Stimulus.png"
width = 30 # in inches
height = 10

#save output
ggsave(plot.fname, width = width, height = height, units="in", limitsize = FALSE)


#### use cowplot to make a grid
####-----------------------------------------------------
p = plot_grid(t.o, ex, st, labels=c('A', 'B', 'C'), label_size = 12)

p # print out the plot so you can see it

# define variables for saving plots
plot.fname = "CC_FX.png"
width = 30 # in inches
height = 10

#save output
p = p + ggsave(plot.fname, width = width, height = height, units="in", limitsize = FALSE)
