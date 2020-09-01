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


# load data and wrangle into tidy form (see https://r4ds.had.co.nz/tidy-data.html), plus relabel to make
# labels a little simpler
dat = read.csv("../total_of_313_subs_AB_task_trial_level_data.csv", header=TRUE)



##### ------------------
task.order.T2gT1 <- dat %>% group_by(Subj.No, Task.Order, Trial.Type) %>%
                            summarise(acc=mean(T2T1.Accuracy))
task.order.T2gT1$Trial.Type <- factor(task.order.T2gT1$Trial.Type)
task.order.T2gT1$Task.Order <- factor(task.order.T2gT1$Task.Order)
task.order.T2gT1$Subj.No <- factor(task.order.T2gT1$Subj.No)

task.order.T2gT1 %>% ggplot(aes(x=Trial.Type, y=acc, fill=Task.Order)) +
                     geom_boxplot() +
                     facet_wrap(~Task.Order)
