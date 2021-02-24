### written by K. Garner, April 2020
### edited by Z. Nott, September 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

### custom functions for running statistical analyses and plotting 

# ----------------------------------------------------------------------------------------------------
# select a data subset
# ----------------------------------------------------------------------------------------------------
sample.N <- function(subs, N, k, replace){
  # get k x N subject numbers and assign to a dataframe for future filtering  
  this.sample <- data.frame(sub=unlist(lapply(N, sample, x=subs, replace=replace)),
                            Nsz=unlist(lapply(N,  function(x) rep(x, each=x))),
                            perm=k)
  this.sample
}

# ----------------------------------------------------------------------------------------------------
###### t-test functions (paired samples)
###### ----------------------------------------------------------------------------

get.ps.t.test <- function(data, iv, dv, x){
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

get.cohens.d <- function(data, iv, dv, x){
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

run.t.test.sim <- function(data, iv="Block.No.Names", dv="RT", x="Random Block", subs){
  
  out = data.frame( p    = get.ps.t.test(data, iv, dv, x),
                    d    = get.cohens.d(data, iv, dv, x))
  out
}

###### -----------------------------------------------------------------------------------------------
###### t-test functions (one sample t test, used for VSL)
###### -----------------------------------------------------------------------------------------------

get.os.t.test <- function(sum.data, dv){
  # run one sample t.test 
  # return the p value
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = iv
  t.dat = sum.data[eval(dv)]
#  t.idx = data[eval(iv)] == x
#  t = t.test(t.dat[t.idx == TRUE], alt = "greater", mu = 0.5)
  t = t.test(t.dat, alt = "greater", mu = 0.5)
  t$p.value
}

get.os.cohens.d <- function(sum.data, dv){
  # get Cohen's d measure for one sample t test
  # sum.data = dataframe for testing
  # dv = dv of interest

  d.dat = sum.data[eval(dv)]
  meanH0 = 0.5
  sd = sd( d.dat$acc )
  d = (mean(d.dat$acc) - meanH0) / sd
  d
}

run.os.t.test.sim <- function(data, dv = "acc"){
  
  sum.data = data %>% group_by(sub) %>%
                      summarise(acc=mean(Response==Target.Order))
  out = data.frame( p    = get.os.t.test(sum.data, dv),
                    d    = get.os.cohens.d(sum.data, dv))
  out
}

###### -----------------------------------------------------------------------------------------------
###### prevalence statistics functions
###### -----------------------------------------------------------------------------------------------

# run.mont.frst.lvl <- function(data, N){
#   # data = 1 participants VSL data! N = number of montecarlo simulations
#   # https://arxiv.org/pdf/1512.00810.pdf
#   # see Algorithm section
#   # As 24! is in the millions, going to do a monte carlo sampling for the first level permutation
#   sub.data <- data.frame(sub=rep(data$sub[1], times=N),
#                          acc=NA)
#   sub.data$acc[1] = with(data, mean(Target.Order==Response))
#   sub.data$acc[2:N]=replicate(N-1, with(data, mean(sample(Target.Order)==Response)))
#   sub.data
# }

run.sing.mont.frst.lvl <- function(data){
  data %>% group_by(sub) %>%
    summarise(acc = mean(sample(Target.Order) == Response))
}

run.mont.frst.lvl.over.subs <- function(data,NfL){
  # feed in all VSL data and the number of monte carlo perms to run (N)
  # will apply the run.mont.frst.lvl over each subject and return a dataframe
  # note: the first permutation is the preserved orderings, as they occurred in the experiment
  nsubs = data$Nsz[1]
  idx <- gen.flvl.perms(NfL, nsubs, max(data$Trial.No))
  names(data)[names(data) == "Trial.No"] = "trial"
  
  perms <- rbind(data %>% group_by(sub) %>%
                          summarise(acc=mean(Response==Target.Order)) %>%
                          mutate(p=1),
                 inner_join(idx, data, by=c("sub", "trial")) %>% group_by(sub, p) %>%
                          summarise(acc=mean(Response==Target.Order)))
  perms %>% arrange(sub)
}

gen.flvl.perms <- function(NfL, nsubs, ntrials, neut=FALSE){
  if (neut == FALSE){
    idx <- data.frame(trial = unlist(replicate(nsubs*(NfL-1), sample(ntrials, replace=FALSE), simplify=FALSE)),
                      p = rep(2:NfL, each=ntrials*nsubs),
                      sub = rep(1:nsubs, each=ntrials, times=NfL-1))
  } else {
    idx <- data.frame(trial = rep(c(1:ntrials), times=nsubs*(NfL-1)),
                      p = rep(2:NfL, each=ntrials*nsubs),
                      sub = rep(1:nsubs, each=ntrials, times=NfL-1))
  }
  idx
}

gen.samps <- function(k, nsubs){
  data <- data.frame(shuffle=rep(2:k, nsubs),
                     sub=rep(1:nsubs, each=k-1),
                     p=unlist(lapply(1:nsubs, function (x) sample(2:k, k-1, replace=T))))
  data
}

get.perm.mins <- function(slvl.idx, flvl.idx, flvl.neut, samp.data){   
  inner_join(slvl.idx, cbind(inner_join(flvl.idx, samp.data, by=c("sub","trial")) %>% select(-Response), inner_join(flvl.neut, samp.data, by=c("sub", "trial")) %>% select(Response)), by=c("sub", "p")) %>%
    group_by(shuffle, sub) %>%
    summarise(acc=mean(Target.Order == Response)) %>%
    group_by(shuffle) %>%
    summarise(min.acc=min(acc)) %>%
    select(min.acc)
}

get.flvl.perms <- function(data, NfL){
  nsubs = max(data$sub)
  data.frame(acc=unlist(lapply(unique(data$sub), 
                               function(x, y=NfL-1) replicate(y, with(data[data$sub == x, ], mean(sample(Target.Order, replace=FALSE) == Response))))),
             sub=rep(1:nsubs, each=NfL-1),
             p = rep(2:NfL, times=nsubs))
}

prev.test <- function(samp.data, alpha, k, NfL){
  # samp.data = sample of data
  # k = the number of 2nd level perms
  # NfL = the number of first level perms
  # alpha = criteria for significance
  
  # first, run first level permutations on input dataset
  # assign a unique subject number to each data entry and
  # get the data from the first level perms
  nsubs <- length(samp.data$Trial.No)/max(samp.data$Trial.No)
  samp.data$sub <- rep(1:nsubs, each=max(samp.data$Trial.No))
  names(samp.data)[names(samp.data) == "Trial.No"] = "trial"
#  flvl.perms <- run.mont.frst.lvl.over.subs(samp.data, NfL) #P1 in 10.1016/j.neuroimage.2016.07.040 (label shuffle)
  # compute flvl idx
  ntrials = max(samp.data$trial)
  flvl.perms <- get.flvl.perms(samp.data, NfL)
  # Now generate a neutral idx as long and as appropriate as flvl.idx, for a subsequent cbind, prior to 
  # summarising and joining to slvl idx, summarising and then taking minimum stat.
  # Now the data is prepared generate indexing for minimum stat
  slvl.idx = rbind(gen.samps(k, nsubs)) %>% arrange(shuffle, sub)
  # computes prevalence statistic, given a set of second level permutations (and original scores)
  # Based on: https://github.com/allefeld/prevalence-permutation/blob/master/prevalenceCore.m - lines 160-168, also
  # k = the number second level permutation you want to extract from the data
  # first select the minimum statistic from the neutral permutation
  # sort out this one 
  neut_m <- min(unlist(lapply(unique(samp.data$sub), function(x) with(samp.data[samp.data$sub == x, ], mean(Response==Target.Order)))))
  # now compute the probability of the minimum value (equation 24 of 10.1016/j.neuroimage.2016.07.040)
 # perm_mins <- t(do.call(rbind, lapply(c(1:max(data$k)), function(x) min(data$acc[data$k == x]))))
  # THE BELOW ONLY WORKS BECAUSE THE DATA ARE IN PROPORTION CORRECT. WOULD NOT WORK FOR PERCENTAGES!
  puGN <- sum(unlist(lapply(unique(slvl.idx$shuffle), function(x) min(inner_join(slvl.idx[slvl.idx$shuffle==x, ], flvl.perms, by=c("sub", "p"))))) >= neut_m)/NfL
  # probability uncorrected of global null (puGN)
  #puGN <- sum(perm.mins$min.acc >= neut_m)/NfL # this is the uncorrected p value for the global null hypothesis that a == a0
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
  
  results <- data.frame(d=gamma_zero,
                        p = puGN)
  results
}

run.prev.test <- function(data, alpha=.05, k=1000, Np=1000){
  # data = the input data for that N and sample
  # N = the desired sample size
  # subs = the list of subjects
  # alpha = the alpha level against which to assess significance
  # k = the number of 2nd level perms
  # Np = the number of 1st level perms
  results <- prev.test(data, alpha, k, Np)
  results
}

# ----------------------------------------------------------------------------------------------------
###### LME and sim functions for SRT data
#### -------------------------------------------------------------------------------------------------
run.lme.4.srt <- function(data, dv = "RT"){
  # run a linear mixed effects analysis
  # going to take the interaction of block x sequence/random as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  
  names(data) <- c("sub", "Nsz", "perm", "trialtype", "RT")  
  
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

# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for CC data
#### -------------------------------------------------------------------------------------------------

get.ps.CC <- function(data, dv = "RT"){
  # run aov on the contextual cueing data
  # return the p value, and the conversion of peta to d
  # data = dataframe for testing
  # dv = name of dv (typically RT)
  
  # NOTE: This also works for SRT data
  
  names(data) <- c("sub", "Nsz", "perm", "block", "trialtype", "RT")
  #an <- aov(eval(parse(text=dv)) ~ (block*trialtype)+Error(sub/(block*trialtype)), data = data) # not worried about using type 1 sum of squares because the data are balanced, see https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
  # comparisons between ours and ss3 outputs revealed the above comment to not be true, therefore using Anova from car package
  an <- Anova(lm(RT ~ block * trialtype, data=data), type=3)
  p <- an$`Pr(>F)`[which(rownames(an) == "block:trialtype")]
  out = list()
  out$p <- p
  # compute partial eta squared: https://www.frontiersin.org/articles/10.3389/fpsyg.2013.00863/full equation 12
  # and convert to d (see https://www.journalofcognition.org/articles/10.5334/joc.10/)
  # peta <- c(summary(an$`Within`)[[1]]["Sum Sq"][2,1]/sum(summary(an$`Within`)[[1]]["Sum Sq"]),
  #           summary(an$`Within`)[[1]]["Sum Sq"][3,1]/sum(summary(an$`Within`)[[1]]["Sum Sq"]))
  peta <-  an$`Sum Sq`[which(rownames(an) == "block:trialtype")]/sum(an$`Sum Sq`[which(rownames(an) == "block:trialtype")],
                                                                     an$`Sum Sq`[which(rownames(an) == "Residuals")])
  out$d <- sqrt((4*peta)/(1-peta))
  out
}


run.lme.4.cc <- function(data, dv = "RT", fx="int"){
  # run a linear mixed effects analysis
  # going to take the interaction of block x repeat/novel as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  # fx = do you want to test for the interaction (int), or the main effect (me)?
  
  names(data) <- c("sub", "Nsz", "perm", "block", "trialtype", "RT")  
  
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


# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for SD data
#### -------------------------------------------------------------------------------------------------
get.ps.SD <- function(data, dv){
  # run aov on the contextual cueing data
  # return the p value, and the conversion of peta to d
  # data = dataframe for testing
  # dv = name of dv (typically RT)
  
  # NOTE: This also works for SRT data
  
  names(data) <- c("sub", "Nsz", "perm", "task", "trialtype", "RT")
  #an <- aov(eval(parse(text=dv)) ~ (task*trialtype)+Error(sub/(task*trialtype)), data = data) # not worried about using type 1 sum of squares because the data are balanced, see https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
  # happy with the defaults
  an <- Anova(lm(RT ~ task * trialtype, data=data), type=3)
  p <- an$`Pr(>F)`[which(rownames(an) == "trialtype")]
  #p <- summary(an$`Within`)[[1]][["Pr(>F)"]][2] # we care about main effect of trial type
  out = list()
  out$p <- p
  # compute partial eta squared and convert to d (see https://www.journalofcognition.org/articles/10.5334/joc.10/)
  peta <-  an$`Sum Sq`[which(rownames(an) == "trialtype")]/sum(an$`Sum Sq`[which(rownames(an) == "trialtype")],
                                                               an$`Sum Sq`[which(rownames(an) == "Residuals")])
  out$d <- sqrt((4*peta)/(1-peta))
  out
}

run.lme.4.SD <- function(data, dv = "RT", efx = "sub"){
  # run a linear mixed effects analysis
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  names(data) <- c("sub", "Nsz", "perm", "task", "trialtype", "RT")
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

# ----------------------------------------------------------------------------------------------------
###### aov and LME functions for AB data
#### -------------------------------------------------------------------------------------------------

get.ps.aov.AB <- function(data, dv="T2gT1"){
  # run aov on the dv with lag as the iv 
  # return the p value
  # data = dataframe for testing - has 4 columns - sub, lag, T1, T2gT1
  # dv = name of dv (T1 or T2gT1)
  names(data) <- c("sub", "Nsz", "perm", "lag", "T1", "T2gT1")
#  an <- aov(eval(parse(text=dv)) ~ (lag)+Error(sub/(lag)), data = data) # not worried about using type 1 sum of squares because the data are balanced, see https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
  an <- Anova(lm(eval(parse(text=dv)) ~ lag, data=data, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
  p <- an$`Pr(>F)`[which(rownames(an) == "lag")]
#  p <- summary(an$`Within`)[[1]][["Pr(>F)"]][1]
  out = list()
  out$p <- p
  # compute partial eta squared and convert to d (see https://www.journalofcognition.org/articles/10.5334/joc.10/)
  # is actually eta squared in this case
  #peta <- summary(an$`Within`)[[1]]["Sum Sq"][1,1]/sum(summary(an$`Within`)[[1]]["Sum Sq"])
  peta <-  an$`Sum Sq`[which(rownames(an) == "lag")]/sum(an$`Sum Sq`[which(rownames(an) == "lag")],
                                                         an$`Sum Sq`[which(rownames(an) == "Residuals")])
  out$d <- sqrt((4*peta)/(1-peta))
  out
}

run.lme.4.AB <- function(data, dv="T2gT1"){
  # run a linear mixed effects analysis
  # going to take the difference between lag 2 and lag 7 as the effect from 
  # which to calculate the recommended d value, 
  # for comparison with the ffx analysis: see https://www.journalofcognition.org/articles/10.5334/joc.10/
  # and https://psycnet-apa-org.ezproxy.library.uq.edu.au/fulltext/2014-32656-001.html
  names(data) <- c("sub", "Nsz", "perm", "lag", "T1", "T2gT1")  
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

# ----------------------------------------------------------------------------------------------------
###### functions to run sims ACROSS ALL TASKS
#### -------------------------------------------------------------------------------------------------
run.outer <- function(in.data, subs, N, k, j, cores, fstem,  ffx.f, rfx.f){
  # this function runs the outer permutation loop
  # for 1:j permutations, for a given N
  # arguments:
  # --in.data = dataframe (see notes of get.ps.aov.AB)
  # --subs = the list of all subject numbers
  # --N = the number to sample from subs
  # --k = the total number of inner permutations
  # --j = which outer permutation are we on?
  # --cores = how many cores do you want to use?
  # -- fstem = what do you want the output files from the inner loop to be called?
  # -- ffx.f & rfx.f: are the references to the functions you want to run for this task
  sub.idx = lapply(1:j, function(x) sample.N(subs, N, x, replace=FALSE))
  sub.idx = do.call(rbind, sub.idx)
  # here I mclapply over each 'parent sample' to pass into run.aov.inner
  # using filter is fine, because parents are sampled without replacement
  mclapply(1:j, function(x) run.inner(in.data=filter(in.data, Subj.No %in% sub.idx$sub[sub.idx$perm==x]),
                                       parent.subs=sub.idx$sub[sub.idx$perm==x],
                                       N=N,
                                       k=k,
                                       j=x,
                                       fstem=fstem,
                                       ffx.f=ffx.f,
                                       rfx.f=rfx.f),
           mc.cores=cores)
}

run.inner <- function(in.data, parent.subs, N, k, j, fstem, ffx.f, rfx.f){
  # this function runs the inner permutation loop
  # for 1:k permutations, run.AB.models is run, and the output
  # collated. The results are saved in a binary file in the
  # local directory
  # arguments:
  # --in.data = dataframe (see notes of get.ps.aov.AB)
  # --subs = the list of all subject numbers
  # --N = the number to sample from subs
  # --perm = which permutation we're on
  # --k = the total number of inner permutations
  # --j = which outer permutation are we on?
  # -- fstem: see run.outer
  
  out = replicate(k, run.models(in.data=in.data, subs=parent.subs, N=N, ffx.f=ffx.f, rfx.f=rfx.f), simplify=FALSE)
  out = do.call(rbind, out)
  out$k = rep(1:k, each=2)
  
  # now save the output
  fname=sprintf(fstem, N, j)
  save(out, file=fname)
}

run.models <- function(in.data, subs, N, ffx.f, rfx.f){
  # this function runs 1 simulation 
  # this function will sample the requested N with replacement
  # and then apply the ffx and rfx analysis and output to a dataframe
  # in.data = dataframe 
  # subs = the list of all subject numbers
  # N = the number to sample from subs

  idx = sample.N(subs=subs, N=N, k=1, replace=TRUE) # get the samples for this one permutation
  names(in.data)[names(in.data)=="Subj.No"] = "sub"
  
  ffx.f <- eval(ffx.f)
  tmp.fx = unlist(ffx.f(inner_join(idx, in.data, by="sub")))
  rfx.f <- eval(rfx.f)
  tmp.rfx = unlist(rfx.f(inner_join(idx, in.data, by="sub")))

  if (length(names(tmp.rfx)) > 2) {
    out = data.frame( n    = c(N, N),
                      p    = c(tmp.fx["p"], tmp.rfx["p"]),
                      d    = c(tmp.fx["d"], tmp.rfx["d"]),
                      esub = c(NA, tmp.rfx["esub"]),
                      eRes = c(NA, tmp.rfx["eRes"]),
                      mod = c("ffx", "rfx") )
  } else {
    out = data.frame( n    = c(N, N),
                      p    = c(tmp.fx["p"], tmp.rfx["p"]),
                      d    = c(tmp.fx["d"], tmp.rfx["d"]),
                      mod = c("ffx", "rfx") )    
  }
  out
}

# ----------------------------------------------------------------------------------------------------
# Density generating functions for plotting, and for computing 95 % CI for the FFX/RFX ratio measure
# ----------------------------------------------------------------------------------------------------

dens.across.N <- function(fstem, Ns, j, min, max, spacer, dv, savekey){
  # grab density functions for dv of choice, across all N sizes
  # save to a binary file output
  # INPUTS
  # -- fstem: see add.dens
  # -- Ns: vector of subject Ns for which there are data to grab
  # -- j: number of outer permutations
  # -- min: see add.dens (and beyond)
  # -- max: see add.dens
  # -- spacer: see get.dens
  # -- dv: which dv are we pulling data out for?
  # -- savekey: usually the initials of the paradigm, for saving output
  tmp = lapply(Ns, add.dens, fstem=fstem, j=j, min=min, max=max, dv=dv, spacer=spacer)
  d = do.call(rbind, tmp)
  fname = paste(savekey, dv, "d.RData", sep="_")
  save(d, file=fname)
}


add.dens <- function(fstem, N, j, min, max, spacer, dv){
  # add density functions for the ffx or rfx based values
  # for a given N, across all outer loops. 
  # output is two density functions.
  # INPUTS
  # -- fstem: fstem = filestem to be sprintf'd with N and j
  # -- N: sub sample size
  # -- j: outer loop size
  # -- min: see get.dens
  # -- max: see get.dens
  # -- spacer: see get.dens
  # -- dv: which dv do you want to know about?

  ds <- lapply(1:j, get.dens, fstem=fstem, N=N, min=min, max=max, dv=dv, spacer=spacer)
  ds <- lapply(1:j, function(x) do.call(rbind, ds[[x]])) 
  ds <- Reduce('+', ds)
  out <- data.frame(Nsz = rep(N, each=ncol(ds)*2),
                    mod = rep(c("ffx", "rfx"), each=ncol(ds)),
                    d = c(ds["ffx",], ds["rfx",]),
                    x = seq(min, max, by=abs(max-min)/spacer))
  out
}


get.dens <- function(fstem, N, j, min, max, dv, spacer){
  # this function loads the relevant data based on fstem, and calls
  # gen.dens to compute a density function
  # for the given dv, for the given N, and j, using min and max as 
  # the density range
  # INPUTS
  # -- fstem = filestem to be sprintf'd with N and j
  # -- N = subject sample size
  # -- j = outer loop permutation number
  # -- min value for density range
  # -- max value for density range
  # -- dv = which d for which to compute density
  
  load(sprintf(fstem, N, j))
  ffx = out[out$mod=="ffx",]
  rfx = out[out$mod=="rfx",]
  if (dv == "p"){
    ffx.d <- gen.dens(min, max, spacer=spacer, log(ffx[,eval(dv)]))
    rfx.d <- gen.dens(min, max, spacer=spacer, log(rfx[,eval(dv)]))
  } else if (dv == "d") {
    ffx.d <- gen.dens(min, max, spacer=spacer, ffx[,eval(dv)])
    rfx.d <- gen.dens(min, max, spacer=spacer, rfx[,eval(dv)])
  } else {
    ffx.d <- rep(NA, length(seq(min, max, by=abs(max-min)/spacer)))
    rfx.d <- gen.dens(min, max, spacer=spacer, rfx[,eval(dv)])
  }
  list(ffx = ffx.d, rfx=rfx.d)
}

gen.dens <- function(min, max, spacer = 10000, data){
  # generate an UNORMALISED density function between min and max, 
  # for the vector 'data'
  # INPUTS
  # -- min: the minimum range on the x of the density function
  # -- max: same, but the max
  # -- data: the vector over which you wish to compute density
  x = seq(min, max, by=abs(max-min)/spacer)
  idx <- sapply(data, function(y) which.min(abs(x-y)))
  idx.idx = seq(1, length(idx), by=1)
  coords = t(rbind(idx, idx.idx))
  d = matrix(0, nrow=length(x), ncol=length(idx))
  # http://eamoncaddigan.net/r/programming/2015/10/22/indexing-matrices/
  d[idx + nrow(d) * (idx.idx-1)] = 1
  apply(d, 1, sum)
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
    ggtitle(paste(data$model[1])) +
    theme(axis.title.x = element_text(face = "italic"))
}


plt.ps <- function(data, xlims, rel_min_height){
  # same as plt.fx.sz but for p values.
  data %>% filter(measure=="p") %>%
    ggplot(mapping=aes(x=value, y=n)) + #, fill=stat(x))) +
    geom_density_ridges(scale=2, rel_min_height=rel_min_height, fill=wes_palette("IsleofDogs1")[1], color=wes_palette("IsleofDogs1")[5]) + # 
    theme_ridges() +
    xlab('p') + ylab('N') + theme_cowplot() + xlim(xlims) +
    geom_vline(aes(xintercept=log(.05)), linetype="dashed") +
    guides(fill = FALSE, colour = FALSE) +
    ggtitle(paste(data$model[1])) +
    theme(axis.title.x = element_text(face = "italic"))
}


RmType <- function(string) { # remove 1st label from facet_wrap
  sub("._", "", string)
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
