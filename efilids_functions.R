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
  tmp <- lapply(sample(subs, size=N, replace=TRUE), function(x) data[data$Subjects == x, ])
  data <- do.call(rbind, tmp)
  data
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

run.t.test.sim <- function(data, iv, dv, x, y, subs, N){
  
  data = get.data(data, subs, N)
  out = data.frame( n    = N,
                    p    = get.ps.t.test(data, iv, dv, x, y),
                    d    = get.cohens.d(data, iv, dv, x, y))
  out
}

###### -----------------------------------------------------------------------------------------------
###### t-test functions (one sample t test, used for VSL)
###### -----------------------------------------------------------------------------------------------

get.os.t.test <- function(data, dv){
  # run one sample t.test 
  # return the p value
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = iv
  t.dat = data[eval(dv)]
#  t.idx = data[eval(iv)] == x
#  t = t.test(t.dat[t.idx == TRUE], alt = "greater", mu = 0.5)
  t = t.test(t.dat, alt = "greater", mu = 0.5)
  t$p.value
}

get.os.cohens.d <- function(data, dv){
  # get Cohen's d measure for one sample t test
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = iv
  d.dat = data[eval(dv)]
  #d.idx = data[eval(iv)] == x
  meanH0 = 0.5
  # sd = sd( d.dat[d.idx == TRUE])
  sd = sd( d.dat$acc )
  # d = (mean(d.dat[d.idx == TRUE]) - meanH0) / sd
  d = (mean(d.dat$acc) - meanH0) / sd
  d
}

run.os.t.test.sim <- function(data, dv, subs, N, perm){
  
  data = get.data(data, subs, N)
  out = data.frame( n    = N,
                    p    = get.os.t.test(data, dv),
                    d    = get.os.cohens.d(data, dv))
  out
}

###### -----------------------------------------------------------------------------------------------
###### prevalence statistics functions
###### -----------------------------------------------------------------------------------------------

run.mont.frst.lvl <- function(data, N){
  # data = 1 participants VSL data! N = number of montecarlo simulations
  # https://arxiv.org/pdf/1512.00810.pdf
  # see Algorithm section
  # As 24! is in the millions, going to do a monte carlo sampling for the first level permutation
  sub.data <- data.frame(sub=rep(data$Subj.No[1], times=N),
                         acc=NA)
  sub.data$acc[1] = with(data, mean(Target.Order==Response))
  sub.data$acc[2:N]=replicate(N-1, with(data, mean(sample(Target.Order)==Response)))
  sub.data
}

run.mont.frst.lvl.over.subs <- function(data,N){
  # feed in all VSL data and the number of monte carlo perms to run (N)
  # will apply the run.mont.frst.lvl over each subject and return a dataframe
  # note: the first permutation is the preserved orderings, as they occurred in the experiment
  subs <- unique(data$Subj.No)
  perms <- lapply(subs, function(x) run.mont.frst.lvl(data=data[data$Subj.No==x,],N=N))
  perms <- do.call(rbind, perms)
  perms$p <- seq(1,N, by=1)
  perms
}

get.4.scnd.lvl <- function(data, k, N, sub){
  # this function extracts data from one participant and 1st level permutation (N)
  data=data[data$sub == sub & data$p == sample(c(2:N),1),]
  data$k=k
  data
}

run.scnd.lvl.mc <- function(data, k, N){
  # this will generate a set of second level permutations across all subjects
  # as defined by https://github.com/allefeld/prevalence-permutation/blob/master/prevalenceCore.m (largely lines 136-157)
  # procedure goes as:
  # 1. for neutral perm, select actual subject data
  # 2. for second level perms, then randomly select a permutation from each subject
  # 3. store the results
  
  # data = the output data from run.mont.frst.lvl.over.subs
  # k = the number of second level permutations
  # N = the number of 1st level permutations
  nsubs = length(unique(data$sub))
  # for each permutation I take a random selection of non-neutral results (2:N), and select
  # the data to make a dataframe
  tmp=do.call(rbind, do.call(rbind, lapply(unique(data$sub), function(x) lapply(c(2:k), get.4.scnd.lvl, data=data, N=N, sub=x))))
  # then I take the neutral data (original permutation and bind it to the permuted second level data)
  neut.dat <- data[data$p == 1, ]
  neut.dat$k <- 1
  out.data <- rbind(neut.dat, tmp)
  out.data
}

prev.test <- function(data, nsubs, alpha, k, N){
  # data = takes the form of out.data from run.scnd.lvl.mc
  # k = the number of 2nd level perms
  # N = the number of first level perms
  data$Subjects <- rep(1:nsubs, each=N) # give each a unique subject identifier
  data$Subjects <- as.factor(data$Subjects)
  names(data)[names(data) == "Subjects"] = "sub"
  data <- run.scnd.lvl.mc(data, k, N)
  
  # computes prevalence statistic, given a set of second level permutations (and original scores)
  # Based on: https://github.com/allefeld/prevalence-permutation/blob/master/prevalenceCore.m - lines 160-168, also
  # k = the number second level permutation you want to extract from the data
  # first select the minimum statistic from the neutral permutation
  neut_m <- min(data$acc[data$k == 1])
  # now compute the probability of the minimum value (equation 24 of 10.1016/j.neuroimage.2016.07.040)
  perm_mins <- t(do.call(rbind, lapply(c(1:max(data$k)), function(x) min(data$acc[data$k == x]))))
  # probability uncorrected of global null (puGN)
  puGN <- sum(perm_mins >= neut_m)/max(data$k) # this is the uncorrected p value for the global null hypothesis that a == a0
  # the above gives the statement of existance, next step is to evaluate against the prevalence null
  # if this is below alpha, then we would say its significant

  ####### attain upper bound on the gamma null that can be rejected (see equation 20 of 10.1016/j.neuroimage.2016.07.040)
  gamma_zero = (alpha^(1/nsubs) - puGN^(1/nsubs)) / (1 - puGN^(1/nsubs))
  if (puGN > alpha) gamma_zero = 0 # not significant so not defined with a prevalence value
  
  ####### the below completes step 5a in the algorithm section of 10.1016/j.neuroimage.2016.07.040
  # # first define prevalence nulls (line 196 or equation 19)
  # null_gammas <- seq(0.01, 1, by = .01)
  # # probability uncorrected of prevalence null
  # puPN <- ((1 - null_gammas) * puGN ^ (1/nsubs) + null_gammas) ^ nsubs
  # sigMN = round(puPN,2) <= alpha
  # # return info in a dataframe
  # if (length(puPN[sigMN]) < 1){
  #   gamma_zero <- 0
  # } else {
  #   gamma_zero <- max(null_gammas[sigMN])
  # }
  
  results <- data.frame(gamma=gamma_zero,
                        p = puGN)
  results
}

run.prev.test <- function(data, N, subs, alpha, k, Np){
  # data = the data output by run.mont.frst.lvl.over.subs
  # N = the desired sample size
  # subs = the list of subjects
  # alpha = the alpha level against which to assess significance
  # k = the number of 2nd level perms
  # Np = the number of 1st level perms
  tmp <- get.data(data, subs, N) #flvl perms dat
  results <- prev.test(tmp, N, alpha, k, Np)
  results$n <- N
  results
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
  
  names(data) <- c("sub", "trialtype", "RT")  
  
  # FOR THE SRT TASK, ALSO FOUND THAT A BOUNDARY (SINGULAR) FIT WHEN INCLUDING RFX OF EXPERIMENTOR AND TASK ORDER
  mod <- lmer(eval(parse(text=dv)) ~ trialtype + (1|sub),
              data=data, REML=FALSE)
  null <- lmer(eval(parse(text=dv)) ~ (1|sub),
               data=data, REML=FALSE)
  d <- abs(summary(mod)$coefficients["trialtypeSequence Block","Estimate"])/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  p <-anova(mod, null)$`Pr(>Chisq)`[2]
  out = list()
  out$p = p
  out$d = d
  df = as.data.frame(VarCorr(mod))
  out$esub = df$sdcor[df$grp=="sub"]^2
#  out$eexp = df$sdcor[df$grp=="exp"]^2
#  out$etask = df$sdcor[df$grp=="task.order"]^2  
  out$eRes = df$sdcor[df$grp=="Residual"]^2  
  out
}

run.SRT.sim <- function(data, subs, N, dv, fx, perm){
  # this function runs the sims for contextual cueing
  # data = conforms to the requirements for get.ps.aov.CC or run.lme.4.cc

  if (fx == "ffx"){
    tmp = run.t.test.sim(data, "Block.No.Names", dv, "Random Block", "Sequence Block", subs, N)
    tmp$esub = NA
#    tmp$etask = NA
#    tmp$eexp = NA
    tmp$eRes = NA
  } else if (fx == "rfx"){
    data = get.data(data, subs, N)
    tmp = run.lme.4.srt(data, dv)
  }
  out = data.frame( n    = N,
                    p    = tmp$p,
                    d    = tmp$d,
                    esub = tmp$esub,
#                    etask = tmp$etask,
#                    eexp = tmp$eexp,
                    eRes = tmp$eRes)
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
  p <- c(summary(an$`Within`)[[1]][["Pr(>F)"]][2], summary(an$`Within`)[[1]][["Pr(>F)"]][3])
  out = list()
  out$p <- p
  # compute partial eta squared and convert to d (see https://www.journalofcognition.org/articles/10.5334/joc.10/)
  peta <- c(summary(an$`Within`)[[1]]["Sum Sq"][2,1]/sum(summary(an$`Within`)[[1]]["Sum Sq"]),
            summary(an$`Within`)[[1]]["Sum Sq"][3,1]/sum(summary(an$`Within`)[[1]]["Sum Sq"]))
  
  out$d <- sqrt((4*peta)/(1-peta))
  out
}

# run.lme.4.cc <- function(data, dv){
#   # run a linear mixed effects analysis
#   # going to take the interaction of block x repeat/novel as the effect from 
#   # which to calculate the recommended d value, 
#   # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
#   # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
#   names(data) <- c("sub", "task.order", "exp", "block", "trialtype", "RT")  
#   mod <- lmer(eval(parse(text=dv)) ~ block*trialtype + (1|sub),
#               data=data, REML=FALSE)
#   null <- lmer(eval(parse(text=dv)) ~ block + (1|sub),
#                data=data, REML=FALSE)
#   d <- abs(summary(mod)$coefficients["block:trialtypeRepeated","Estimate"])/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
#   p <-anova(mod, null)$`Pr(>Chisq)`[2]
#   out = list()
# 
#   df = as.data.frame(VarCorr(mod))
#   out$esub = df$sdcor[df$grp=="sub"]^2
#   # out$eexp = df$sdcor[df$grp=="exp"]^2
#   # out$etask = df$sdcor[df$grp=="task.order"]^2  
#   out$eRes = df$sdcor[df$grp=="Residual"]^2  
#   
#   out$p = p
#   out$d = d
#   out
# }

run.lme.4.cc <- function(data, dv, fx){
  # run a linear mixed effects analysis
  # going to take the interaction of block x repeat/novel as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  # fx = do you want to test for the interaction (int), or the main effect (me)?
  
  names(data) <- c("sub", "task.order", "exp", "block", "trialtype", "RT")  
  
  if (fx == "int"){
    mod <- lmer(eval(parse(text=dv)) ~ block + trialtype + block:trialtype + (1|sub),
                data=data, REML=FALSE)
    null <- lmer(eval(parse(text=dv)) ~ block + trialtype + (1|sub),
                 data=data, REML=FALSE)
    d <- abs(summary(mod)$coefficients["block:trialtypeRepeated","Estimate"])/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  } else if (fx == "me"){
    mod <- lmer(eval(parse(text=dv)) ~ block + trialtype + (1|sub),
                data=data, REML=FALSE)
    null <- lmer(eval(parse(text=dv)) ~ block + (1|sub),
                 data=data, REML=FALSE)
    d <- abs(summary(mod)$coefficients["trialtypeRepeated","Estimate"])/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  }
  
  p <- anova(mod, null)$`Pr(>Chisq)`[2]
  out = list()
  
  df = as.data.frame(VarCorr(mod))
  out$esub = df$sdcor[df$grp=="sub"]^2
  # out$eexp = df$sdcor[df$grp=="exp"]^2
  # out$etask = df$sdcor[df$grp=="task.order"]^2  
  out$eRes = df$sdcor[df$grp=="Residual"]^2  
  
  out$p = p
  out$d = d
  out
}


run.aov.CC.sim <- function(data, subs, N, dv, fx, efx){
# this function runs the sims for contextual cueing
# data = conforms to the requirements for get.ps.aov.CC or run.lme.4.cc
# efx = the effect you want to test for in the CC data 'me' or 'int' 
  data = get.data(data, subs, N)
  if (fx == "ffx"){

        tmp = get.ps.aov.CC.SRT(data, dv)
    tmp$esub = NA
    # tmp$etask = NA
    # tmp$eexp = NA
    tmp$eRes = NA
    tmp$fx = c("me", "int")
    
  } else if (fx == "rfx"){
    
    tmp = run.lme.4.cc(data, dv, efx)
    tmp$fx = paste(efx)
  }
  out = data.frame( n    = N,
                    p    = tmp$p,
                    d    = tmp$d,
                    esub = tmp$esub,
                    # etask = tmp$etask,
                    # eexp = tmp$eexp,
                    eRes = tmp$eRes,
                    fx = tmp$fx)
  out
}


# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for SD data
#### -------------------------------------------------------------------------------------------------
get.ps.aov.SD <- function(data, dv){
  # run aov on the contextual cueing data
  # return the p value, and the conversion of peta to d
  # data = dataframe for testing
  # dv = name of dv (typically RT)
  
  # NOTE: This also works for SRT data
  
  names(data) <- c("sub", "task", "trialtype", "RT")
  an <- aov(eval(parse(text=dv)) ~ (task*trialtype)+Error(sub/(task*trialtype)), data = data) # not worried about using type 1 sum of squares because the data are balanced, see https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
  p <- summary(an$`Within`)[[1]][["Pr(>F)"]][2] # we care about main effect of trial type
  out = list()
  out$p <- p
  # compute partial eta squared and convert to d (see https://www.journalofcognition.org/articles/10.5334/joc.10/)
  peta <- summary(an$`Within`)[[1]]["Sum Sq"][2,1]/sum(summary(an$`Within`)[[1]]["Sum Sq"])
  out$d <- sqrt((4*peta)/(1-peta))
  out
}

run.lme.4.SD <- function(data, dv, efx){
  # run a linear mixed effects analysis
  # going to take the difference between lag 2 and lag 7 as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  if (sum(names(data) == "task.stim")==0 ){
    names(data) <- c("sub", "task", "trialtype", "RT")  
  } else {
    names(data) <- c("sub", "task", "trialtype", "stim", "RT")  
  }
  # NOTE: when I had the following RFX structure (1|sub) + (1|task.order) + (1|exp) the fit was singular, suggesting
  # the RFX structure is too complex for the data - https://stats.stackexchange.com/questions/378939/dealing-with-singular-fit-in-mixed-models
  # therefore am testing the removal of each element of the rfx structure, starting with task.order.
  # Dropping task.order (with N=23) was far superior to dropping experimenter:
  # Formula: eval(parse(text = dv)) ~ trialtype * task + (1 | sub) + (1 |      exp)
  # Data: tmp
  # AIC        BIC     logLik   deviance   df.resid 
  # -1675.0323 -1639.1048   844.5161 -1689.0323       1245 
  
  # Formula: eval(parse(text = dv)) ~ trialtype * task + (1 | sub) + (1 |      task.order)
  # Data: data
  # AIC       BIC    logLik  deviance  df.resid 
  # -132.5054 -114.8529   73.2527 -146.5054        85 
  # FITS still singular with exp, so dropping that also. The remaining model only has a singular fit a small number of times
  if (efx == "sub"){
    mod <- lmer(eval(parse(text=dv)) ~ trialtype*task + (1|sub),
                data=data, REML=FALSE)
    null <- lmer(eval(parse(text=dv)) ~ task + (1|sub),
                 data=data, REML=FALSE)
   
  } else if (efx == "stim"){
    mod <- lmer(eval(parse(text=dv)) ~ trialtype*task + (1|sub) + (1|stim),
                data=data, REML=FALSE)
    null <- lmer(eval(parse(text=dv)) ~ task + (1|sub) + (1|stim),
                 data=data, REML=FALSE)
  }
  
  d <- abs(summary(mod)$coefficients["trialtypesingle","Estimate"])/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  p <- anova(mod, null)$`Pr(>Chisq)`[2]
  out = list()
  out$p = p
  out$d = d
  # below is hard coded and clunky - beware! largely done this way because of past decisions!
  df = as.data.frame(VarCorr(mod))
  out$esub = df$sdcor[df$grp=="sub"]^2
  out$eRes = df$sdcor[df$grp=="Residual"]^2
  if (efx == "stim") out$estim = df$sdcor[df$grp=="stim"]^2
  #list(out, as.data.frame(VarCorr(mod)))
  out
}

run.SD.sim <- function(data, subs, N, dv, fx, efx){
  # this function runs the sims for SD
  # data = conforms to the requirements for get.ps.aov.SD or run.lme.4.SD
  data = get.data(data, subs, N)
  if (fx == "ffx"){
    tmp = get.ps.aov.SD(data, dv)
    tmp$esub = NA
    tmp$eexp = NA
    tmp$eRes = NA
  } else if (fx == "rfx"){
    tmp = run.lme.4.SD(data, dv, efx)
  }
  out = data.frame( n    = N,
                    p    = tmp$p,
                    d    = tmp$d,
                    esub = tmp$esub,
                    eRes = tmp$eRes,
                    fx   = efx)
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
  mod <- lmer(eval(parse(text=dv)) ~ lag + (1|sub),
              data=data, REML=FALSE)
  null <- lmer(eval(parse(text=dv)) ~ (1|sub),
               data=data, REML=FALSE)
  d <- summary(mod)$coefficients["laglag_7","Estimate"]/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  p <-anova(mod, null)$`Pr(>Chisq)`[2]
  out = list()
  
  df = as.data.frame(VarCorr(mod))
  out$esub = df$sdcor[df$grp=="sub"]^2
  # out$eexp = df$sdcor[df$grp=="exp"]^2
  # out$etask = df$sdcor[df$grp=="task.order"]^2  
  out$eRes = df$sdcor[df$grp=="Residual"]^2  
  out$p = p
  out$d = d
  
  out
}

run.stim.lme.4.aov <- function(data, dv){
  # run a linear mixed effects analysis
  # going to take the difference between lag 2 and lag 7 as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  
  # notes on model comparison and selection
  # 1. mod: eval(parse(text = dv)) ~ lag + (1 | sub)
  # T1int: eval(parse(text = dv)) ~ lag + (1 | sub) + (1 | T1stim)
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod      6 30714 30764 -15351    30702                         
  # T1int    7 30496 30554 -15241    30482 219.88  1  < 2.2e-16 ***
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 2/ Models:
  # T1int: eval(parse(text = dv)) ~ lag + (1 | sub) + (1 | T1stim)
  # T2int: eval(parse(text = dv)) ~ lag + (1 | sub) + (1 | T1stim) + (1 | 
  #                                                                     T2int:     T2stim)
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # T1int    7 30496 30554 -15241    30482                         
  # T2int    8 30006 30072 -14995    29990 492.38  1  < 2.2e-16 ***
  # any models with random slopes showed a boundary (singular fit)
  
  names(data) <- c("sub", "lag", "T1stim", "T1", "T2gT1")  

  mod <- lmer(eval(parse(text=dv)) ~ lag + (1|sub) + (1|T1stim),
                data=data, REML=FALSE)
  null <- lmer(eval(parse(text=dv)) ~ (1|sub) + (1|T1stim),
               data=data, REML=FALSE)
  d <- summary(mod)$coefficients["laglag_7","Estimate"]/sqrt(sum(as.data.frame(VarCorr(mod))$sdcor^2)) # get the variance of the random effects
  p <-anova(mod, null)$`Pr(>Chisq)`[2]
  out = list()
  
  df = as.data.frame(VarCorr(mod))
  out$esub = df$sdcor[df$grp=="sub"]^2
  out$T1stim = df$sdcor[df$grp=="T1stim"]^2
  out$eRes = df$sdcor[df$grp=="Residual"]^2  
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
  if (fx != "stimrfx"){
  if (fx == "ffx"){
    tmp = get.ps.aov.AB(data, dv)
    tmp$esub = NA
    # tmp$etask = NA
    # tmp$eexp = NA
    tmp$eRes = NA
  } else if (fx == "rfx"){
    tmp = run.lme.4.aov(data, dv)
  }
  out = data.frame( n    = N,
                    p    = tmp$p,
                    d    = tmp$d,
                    esub = tmp$esub,
                    # etask = tmp$etask,
                    # eexp = tmp$eexp,
                    eRes = tmp$eRes)
  } else {
    tmp = run.stim.lme.4.aov(data, dv)
    out = data.frame(n = N,
                     p = tmp$p,
                     d = tmp$d,
                     eSub = tmp$esub,
                     eT1stim = tmp$T1stim,
                     eRes = tmp$eRes)
  }
  out
}


# ----------------------------------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------------------------------

plt.fx.sz <- function(data, ylims){
  # plot effect size, given dataframe of 'n', 'measure', and 'value'
  data %>% filter(measure=="d") %>%
    ggplot(mapping=aes(x=value, y=n)) + #, fill=stat(x))) +
    geom_density_ridges(scale=2, rel_min_height=.01, fill=wes_palette("IsleofDogs1")[1], color=wes_palette("IsleofDogs1")[5]) +
#    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
    theme_ridges() +
#    scale_fill_viridis_c(name = "value", option = "C") +
    xlab('d') + ylab('N') + theme_cowplot() + xlim(ylims) +
 #   scale_y_discrete(breaks = seq(23, 303, by = 20), labels=as.character(seq(23, 303, by = 20))) +
    guides(fill = FALSE, colour = FALSE) +
    ggtitle(paste(data$model[1], data$fx[1], sep=" ")) +
    theme(axis.title.x = element_text(face = "italic"))
}


plt.ps <- function(data, xlims){
  # same as plt.fx.sz but for p values.
  data %>% filter(measure=="p") %>%
    ggplot(mapping=aes(x=value, y=n)) + #, fill=stat(x))) +
    geom_density_ridges(scale=2, rel_min_height=.01, fill=wes_palette("IsleofDogs1")[1], color=wes_palette("IsleofDogs1")[4]) +
    theme_ridges() +
    xlab('p') + ylab('N') + theme_cowplot() + 
    xlim(xlims) +
    geom_vline(aes(xintercept=.05), linetype="dashed") +
    guides(fill = FALSE, colour = FALSE) +
    ggtitle(paste(data$model[1], data$fx[1], sep=" ")) +
    theme(axis.title.x = element_text(face = "italic"))
}


plt.rfx <- function(data, xlims){
  # same as plt.fx.sz but for p values.
  data %>% pivot_longer(names(data)[!names(data) %in% c("n","model")], names_to = "rfx", values_to="var") %>%
    drop_na() %>%
    ggplot(mapping=aes(x=var, y=n)) + #, fill=stat(x))) +
    geom_density_ridges(scale=2, rel_min_height=.01, fill=wes_palette("IsleofDogs1")[5], color=wes_palette("IsleofDogs1")[4]) +
    #    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
    theme_ridges() +
#    scale_fill_viridis_c(name = "value", option = "C") +
    xlab(expression(sigma)) + ylab('N') + theme_cowplot() + xlim(xlims) + 
#    scale_y_discrete(breaks = seq(23, 303, by = 20), labels=as.character(seq(23, 303, by = 20))) +
    facet_wrap(~model*rfx) +
    guides(fill = FALSE, colour = FALSE) +
    theme(axis.title.x = element_text(face = "italic"))
}
