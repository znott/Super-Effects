### written by K. Garner, April 2020
### edited by Z. Nott, September 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

### custom functions for running statistical analyses and plotting 
get.data <- function(data, subs, N){
  # given a dataframe, the list of subject numbers, & N (the number of subs required),
  # return a dataframe with that number of subs returned by replacement
  data[data$Subjects %in% sample(subs, size=N, replace=TRUE), ]
}

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

###### aov functions for contextual cueing
###### ----------------------------------------------------------------------------
get.ps.aov.CC.test <- function(data){
  # run t.test on the dv between the variables x & y 
  # return the p value
  # data = dataframe for testing
  # dv = name of dv

  an <- aov(RT_ms_trimmed ~ (Condition*Epoch)+Error(Subjects/(Condition*Epoch)), data = data) # not worried about using type 1 sum of squares because the data are balanced, see https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
  p <- summary(an$`Subjects:Condition:Epoch`)[[1]][["Pr(>F)"]]
  out = list()
  out$p <- p[!is.na(p)]
  # compute partial eta squared
  out$peta <- summary(an$`Subjects:Condition:Epoch`)[[1]]["Sum Sq"][1,1]/sum(summary(an$`Subjects:Condition:Epoch`)[[1]]["Sum Sq"])
  out
}

run.aov.CC.sim <- function(data, subs, N, perm){
# this function runs the sims using the contextual cueing aov  
  data = get.data(data, subs, N)
  tmp = get.ps.aov.CC.test(data)
  out = data.frame( n    = N,
                    p    = tmp$p,
                    peta  = tmp$peta)
  out
}

###### aov functions for ANOVA
#### ------------------------------------------------------------------
###### aov functions for contextual cueing
###### ----------------------------------------------------------------------------
get.ps.aov.AB.test <- function(data){
  # run t.test on the dv between the variables x & y 
  # return the p value
  # data = dataframe for testing
  # dv = name of dv
  
  an <- aov(acc ~ (lag)+Error(Subjects/(lag)), data = data) # not worried about using type 1 sum of squares because the data are balanced, see https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
  p <- summary(an$`Subjects:lag`)[[1]][["Pr(>F)"]]
  out = list()
  out$p <- p[!is.na(p)]
  # compute partial eta squared
  out$peta <- summary(an$`Subjects:lag`)[[1]]["Sum Sq"][1,1]/sum(summary(an$`Subjects:lag`)[[1]]["Sum Sq"])
  out
}

run.aov.AB.sim <- function(data, subs, N, perm){
  # this function runs the sims using the contextual cueing aov  
  data = get.data(data, subs, N)
  tmp = get.ps.aov.AB.test(data)
  out = data.frame( n    = N,
                    p    = tmp$p,
                    peta  = tmp$peta)
  out
}
