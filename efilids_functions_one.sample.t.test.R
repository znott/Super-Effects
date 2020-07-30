### written by K. Garner, April 2020
### edited by Z. Nott, July 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

### custom functions for running statistical analyses and plotting 
get.data <- function(data, subs, N){
  # given a dataframe, the list of subject numbers, & N (the number of subs required),
  # return a dataframe with that number of subs returned by replacement
  data[data$Subjects %in% sample(subs, size=N, replace=TRUE), ]
}

get.ps.t.test <- function(data, iv, dv, x){
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

get.cohens.d <- function(data, iv, dv, x){
  # get Cohen's d measure for one sample t test
  # data = dataframe for testing
  # iv = name of iv
  # dv = name of dv
  # x = iv
  d.dat = data[eval(dv)]
  d.idx = data[eval(iv)] == x
  meanH0 = 0.5
  sd = sd( d.dat[d.idx == TRUE])
  d = (meanH0 - mean(d.dat[d.idx == TRUE])) / sd
  d
}

run.t.test.sim <- function(data, iv, dv, x, subs, N, perm){
  
  data = get.data(data, subs, N)
  out = data.frame( n    = N,
                    p    = get.ps.t.test(data, iv, dv, x),
                    d    = get.cohens.d(data, iv, dv, x))
  out
}
