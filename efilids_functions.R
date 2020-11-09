### written by K. Garner, April 2020
### edited by Z. Nott, September 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

### custom functions for running statistical analyses and plotting 

# ----------------------------------------------------------------------------------------------------
# select a data subset
# ----------------------------------------------------------------------------------------------------
get.data <- function(data, subs, N){
  # given a dataframe, the list of subject numbers, & N (the number of subs required),
  # return a dataframe with that number of subs returned by replacement
  # **first column of data must be the subject numbers!**
  names(data)[1] = "Subjects"
  data[data$Subjects %in% sample(subs, size=N, replace=TRUE), ]
}

# ----------------------------------------------------------------------------------------------------
###### t-test functions (paired samples)
###### ----------------------------------------------------------------------------

get.ps.t.test <- function(data, iv, dv, x, y){
  # run t.test on the dv between the variables x & y 
  # return the p value
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = first level of iv
  # y = second level of iv
  t.dat = data[eval(dv)]
  t.idx = data[eval(iv)] == x
  t = t.test(t.dat[t.idx == TRUE], t.dat[t.idx==FALSE])
  t$p.value
}

get.cohens.d <- function(data, iv, dv, x, y){
  # get Cohen's d measure for paired samples
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = first level of iv
  # y = second level of iv  
  d.dat = data[eval(dv)]
  d.idx = data[eval(iv)] == x
  sd.pool = sqrt( ( sd( d.dat[d.idx == TRUE]  )^2 + sd(  d.dat[d.idx == FALSE]  )^2 ) / 2 )
  d = (mean(d.dat[d.idx == TRUE]) - mean(d.dat[d.idx == FALSE])) / sd.pool
  d
}

run.t.test.sim <- function(data, iv, dv, x, y, subs, N, perm){
  
  data = get.data(data, subs, N)
  out = data.frame( n    = N,
                    p    = get.ps.t.test(data, iv, dv, x, y),
                    d    = get.cohens.d(data, iv, dv, x, y))
  out
}

###### t-test functions (one sample t test)
###### ----------------------------------------------------------------------------

get.os.t.test <- function(data, iv, dv, x){
  # run one sample t.test 
  # return the p value
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = iv
  t.dat = data[eval(dv)]
  t.idx = data[eval(iv)] == x
  t = t.test(t.dat[t.idx == TRUE], alt = "greater", mu = 0.5)
  t$p.value
}

get.os.cohens.d <- function(data, iv, dv, x){
  # get Cohen's d measure for one sample t test
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = iv
  d.dat = data[eval(dv)]
  d.idx = data[eval(iv)] == x
  meanH0 = 0.5
  sd = sd( d.dat[d.idx == TRUE])
  d = (mean(d.dat[d.idx == TRUE]) - meanH0) / sd
  d
}

run.os.t.test.sim <- function(data, iv, dv, x, subs, N, perm){
  
  data = get.data(data, subs, N)
  out = data.frame( n    = N,
                    p    = get.os.t.test(data, iv, dv, x),
                    d    = get.os.cohens.d(data, iv, dv, x))
  out
}

# ----------------------------------------------------------------------------------------------------
###### LME and sim functions for SRT data
#### -------------------------------------------------------------------------------------------------

# NOTE: the ffx aov model for CC (defined below) also works for SRT

run.lme.4.srt <- function(data, dv){
  # run a linear mixed effects analysis
  # going to take the interaction of block x sequence/random as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  names(data) <- c("sub", "task.order", "exp", "block", "trialtype", "RT")  
  mod <- lmer(eval(parse(text=dv)) ~ block*trialtype + (1|sub) + (1|task.order) + (1|exp),
              data=data, REML=FALSE)
  null <- lmer(eval(parse(text=dv)) ~ block + (1|sub) + (1|task.order) + (1|exp),
               data=data, REML=FALSE)
  d <- abs(summary(mod)$coefficients["block:trialtypeSequence Block","Estimate"])/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  p <-anova(mod, null)$`Pr(>Chisq)`[2]
  out = list()
  out$p = p
  out$d = d
  out
}

run.SRT.sim <- function(data, subs, N, dv, fx){
  # this function runs the sims for contextual cueing
  # data = conforms to the requirements for get.ps.aov.CC or run.lme.4.cc
  data = get.data(data, subs, N)
  if (fx == "ffx"){
    tmp = get.ps.aov.CC.SRT(data, dv)
  } else if (fx == "rfx"){
    tmp = run.lme.4.srt(data, dv)
  }
  out = data.frame( n    = N,
                    p    = tmp$p,
                    d    = tmp$d)
  out
}


# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for CC data
#### -------------------------------------------------------------------------------------------------

get.ps.aov.CC.SRT <- function(data, dv){
  # run aov on the contextual cueing data
  # return the p value, and the conversion of peta to d
  # data = dataframe for testing
  # dv = name of dv (typically RT)
  
  # NOTE: This also works for SRT data
  
  names(data) <- c("sub", "block", "trialtype", "RT")
  an <- aov(eval(parse(text=dv)) ~ (block*trialtype)+Error(sub/(block*trialtype)), data = data) # not worried about using type 1 sum of squares because the data are balanced, see https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
  p <- summary(an$`Within`)[[1]][["Pr(>F)"]][3]
  out = list()
  out$p <- p
  # compute partial eta squared and convert to d (see https://www.journalofcognition.org/articles/10.5334/joc.10/)
  peta <- summary(an$`Within`)[[1]]["Sum Sq"][3,1]/sum(summary(an$`Within`)[[1]]["Sum Sq"])
  out$d <- sqrt((4*peta)/(1-peta))
  out
}

run.lme.4.cc <- function(data, dv){
  # run a linear mixed effects analysis
  # going to take the interaction of block x repeat/novel as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  names(data) <- c("sub", "task.order", "exp", "block", "trialtype", "RT")  
  mod <- lmer(eval(parse(text=dv)) ~ block*trialtype + (1|sub) + (1|task.order) + (1|exp),
              data=data, REML=FALSE)
  null <- lmer(eval(parse(text=dv)) ~ block + (1|sub) + (1|task.order) + (1|exp),
               data=data, REML=FALSE)
  d <- abs(summary(mod)$coefficients["block:trialtypeRepeated","Estimate"])/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  p <-anova(mod, null)$`Pr(>Chisq)`[2]
  out = list()
  out$p = p
  out$d = d
  out
}

run.aov.CC.sim <- function(data, subs, N, dv, fx){
# this function runs the sims for contextual cueing
# data = conforms to the requirements for get.ps.aov.CC or run.lme.4.cc
  data = get.data(data, subs, N)
  if (fx == "ffx"){
    tmp = get.ps.aov.CC.SRT(data, dv)
  } else if (fx == "rfx"){
    tmp = run.lme.4.cc(data, dv)
  }
  out = data.frame( n    = N,
                    p    = tmp$p,
                    d    = tmp$d)
  out
}

# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for AB data
#### -------------------------------------------------------------------------------------------------

get.ps.aov.AB <- function(data, dv){
  # run aov on the dv with lag as the iv 
  # return the p value
  # data = dataframe for testing - has 4 columns - sub, lag, T1, T2gT1
  # dv = name of dv (T1 or T2gT1)
  
  names(data) <- c("sub", "lag", "T1", "T2gT1")
  an <- aov(eval(parse(text=dv)) ~ (lag)+Error(sub/(lag)), data = data) # not worried about using type 1 sum of squares because the data are balanced, see https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
  p <- summary(an$`Within`)[[1]][["Pr(>F)"]][1]
  out = list()
  out$p <- p
  # compute partial eta squared and convert to d (see https://www.journalofcognition.org/articles/10.5334/joc.10/)
  peta <- summary(an$`Within`)[[1]]["Sum Sq"][1,1]/sum(summary(an$`Within`)[[1]]["Sum Sq"])
  out$d <- sqrt((4*peta)/(1-peta))
  out
}

run.lme.4.aov <- function(data, dv){
  # run a linear mixed effects analysis
  # going to take the difference between lag 2 and lag 7 as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  names(data) <- c("sub", "task.order", "exp", "lag", "T1", "T2gT1")  
  mod <- lmer(eval(parse(text=dv)) ~ lag + (1|sub) + (1|task.order) + (1|exp),
              data=data, REML=FALSE)
  null <- lmer(eval(parse(text=dv)) ~ + (1|sub) + (1|task.order) + (1|exp),
               data=data, REML=FALSE)
  d <- summary(mod)$coefficients["laglag_7","Estimate"]/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  p <-anova(mod, null)$`Pr(>Chisq)`[2]
  out = list()
  out$p = p
  out$d = d
  out
}

run.aov.AB.sim <- function(data, subs, N, dv, fx){
  # this function runs the sims using get.ps.aov.AB
  # data = dataframe (see notes of get.ps.aov.AB)
  # subs = the list of all subject numbers
  # N = the number to sample from subs
  # perm = which permutation we're on
  # dv = which DV to analyse (T1 or T2gT1)
  # fx = if ffx run the fixed effects analysis, if rfx, run the rfx analysis
  data = get.data(data, subs, N)
  if (fx == "ffx"){
    tmp = get.ps.aov.AB(data, dv)
  } else if (fx == "rfx"){
    tmp = run.lme.4.aov(data, dv)
  }
  out = data.frame( n    = N,
                    p    = tmp$p,
                    d    = tmp$d)
  out
}


# ----------------------------------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------------------------------

plt.fx.sz <- function(data){
  # plot effect size, given dataframe of 'n', 'measure', and 'value'
  data %>% filter(measure=="d") %>%
    ggplot(mapping=aes(x=n, y=value, fill = n, colour = n)) +
    geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim =
                       TRUE) +
    geom_boxplot(aes(x = as.numeric(n)+0.25, y = value), outlier.shape = NA,
                 alpha = 0.3, width = .1, colour = "BLACK") +
    ylab('d') + xlab('N') + theme_cowplot() + 
    guides(fill = FALSE, colour = FALSE) +
    coord_flip() + ggtitle(data$model[1]) +          
    theme(axis.title.x = element_text(face = "italic"))
}

plt.ps <- function(data){
  # same as plt.fx.sz but for p values.
  data %>% filter(measure == "p") %>%
    ggplot(mapping=aes(x=n, y=value, fill = n, colour = n)) +
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
}