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

# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
dat = read.csv("total_of_313_subs_SRT_task_trial_level_data.csv", header=TRUE)

##### RT.ms: Trial Type x Task Order
##### --------------------------------------------------------
task.order.srt <- dat %>% group_by(Subj.No, Task.Order, Trial.Type, RT.ms) %>% summarise(mean=mean(RT.ms))
task.order.srt$Trial.Type <- factor(task.order.srt$Trial.Type)
task.order.srt$Task.Order <- factor(task.order.srt$Task.Order)
task.order.srt$Subj.No <- factor(task.order.srt$Subj.No)

t.o <- task.order.srt %>% ggplot(aes(x=Task.Order, y=RT.ms, fill=Trial.Type)) + ylim(0,2000) +
  geom_boxplot() +
  facet_wrap(~Trial.Type)
t.o

##### RT_ms: Trial Type x Experimenter
##### --------------------------------------------------------

experimenter.srt <- dat %>% group_by(Subj.No, Experimenter, Trial.Type, RT.ms) %>% summarise(mean=mean(RT.ms))
experimenter.srt$Subj.No <- factor(experimenter.srt$Subj.No)
experimenter.srt$Trial.Type <- factor(experimenter.srt$Trial.Type)

ex <- experimenter.srt %>% ggplot(aes(x=Experimenter, y=RT.ms, fill=Trial.Type,)) + ylim(0,2000) +
  geom_boxplot() +
  facet_wrap(~Trial.Type)
ex

#### use cowplot to make a grid
####-----------------------------------------------------
p = plot_grid(t.o, ex, labels=c('A', 'B'), label_size = 12)

p # print out the plot so you can see it

# define variables for saving plots
plot.fname = "SRT_FX.png"
width = 30 # in inches
height = 10

#save output
p = p + ggsave(plot.fname, width = width, height = height, units="in", limitsize = FALSE)
