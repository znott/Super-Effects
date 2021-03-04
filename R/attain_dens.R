### written by K. Garner, April 2020
### edited by Z. Nott, July 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

# ----------------------------------------------------------------------------------------------------
rm(list=ls())
# ----------------------------------------------------------------------------------------------------
### ATTAIN DENSITIES FOR ANY TASK GIVEN RANGE INPUTS
# ----------------------------------------------------------------------------------------------------
# load packages and source function files
# ----------------------------------------------------------------------------------------------------

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the location of this file
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
# define ref variables
# ----------------------------------------------------------------------------------------------------
sub.Ns = round(exp(seq(log(13), log(313), length.out = 20)))
n.perms =1000# for each sample size, we will repeat our experiment n.perms^2 times 
cores = 30
datpath = "../data/"
rxvnme = "supfxrdat.zip"
task = "AB"

# ----------------------------------------------------------------------------------------------------
# attain densities for each subject N, across all outer samples
# ----------------------------------------------------------------------------------------------------
#dens.across.N(fstem="AB_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=-800, max=0, spacer=1000, dv="p", savekey="AB")
dens.across.N(fstem="_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=10, spacer=500, dv="d", savekey=task, datpath=datpath, rxvnme=rxvnme)
#dens.across.N(fstem="AB_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=0.75, spacer=1000, dv="esub", savekey="AB")
#dens.across.N(fstem="AB_N-%d_parent-%d.RData", Ns=sub.Ns, j=n.perms, min=0, max=0.75, spacer=1000, dv="eRes", savekey="AB")