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
dat = read.csv("../total_of_313_subs_AB_task_trial_level_data.csv", header=TRUE)

##### T2|T1: Trial Type x Task Order
##### --------------------------------------------------------
task.order.T2gT1 <- dat %>% group_by(Subj.No, Task.Order, Trial.Type) %>%
  summarise(acc=mean(T2T1.Accuracy))
task.order.T2gT1$Trial.Type <- factor(task.order.T2gT1$Trial.Type)
task.order.T2gT1$Task.Order <- factor(task.order.T2gT1$Task.Order)
task.order.T2gT1$Subj.No <- factor(task.order.T2gT1$Subj.No)

t.o <- task.order.T2gT1 %>% ggplot(aes(x=Task.Order, y=acc, fill=Trial.Type)) +
  geom_boxplot() +
  facet_wrap(~Trial.Type)

t.o <- task.order.T2gT1 %>% ggplot(aes(x=Trial.Type, y=acc, fill=Trial.Type)) +
  geom_boxplot() +
  facet_wrap(~Task.Order)
t.o

##### T2|T1: Trial Type x Experimenter
##### --------------------------------------------------------

experimenter.T2gT1 <- dat %>% group_by(Subj.No, Experimenter, Trial.Type) %>%
  summarise(acc=mean(T2T1.Accuracy))
experimenter.T2gT1$Subj.No <- factor(experimenter.T2gT1$Subj.No)
experimenter.T2gT1$Trial.Type <- factor(experimenter.T2gT1$Trial.Type)

### plot data
ex <- experimenter.T2gT1 %>% ggplot(aes(x=Experimenter, y=acc, fill=Trial.Type,)) +
  geom_boxplot() +
  facet_wrap(~Trial.Type)
experimenter.T2gT1 %>% ggplot(aes(x=Trial.Type, y=acc, fill=Trial.Type,)) +
  geom_boxplot() +
  facet_wrap(~Experimenter)
ex

# average for each experimenter
experimenter.T2gT1 %>% group_by(Experimenter, Trial.Type) %>%
                       summarise(acc=mean(acc)) %>%
                       ggplot(aes(x=Trial.Type, y=acc, group=Experimenter)) +
                       geom_line()


#### use cowplot to make a grid
####-----------------------------------------------------
p = plot_grid(t.o, ex, labels=c('A', 'B'), label_size = 12)

p # print out the plot so you can see it

# define variables for saving plots
plot.fname = "AB_FX.png"
width = 30 # in inches
height = 10

#save output
p = p + ggsave(plot.fname, width = width, height = height, units="in", limitsize = FALSE)

###---------------------------------------------------------------------
### create barplot to see experimenter frequency for T2|T1 data
dat.ex.T2gT1 <- table(experimenter.T2gT1$Experimenter)
range(dat.ex.T2gT1)
view(dat.ex.T2gT1)
par( mar = c( 2.1, 4.1, 2.1, 2.1) )
barplot(height = dat.ex.T2gT1, ylab = "Frequency", main = "Participants per Experimenter", names.arg = 'Experimenters')
 
