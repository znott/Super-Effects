### written by K. Garner, April 2020
### for the project 'On the detectability of effects in executive function and implicit learning tasks'
### Garner, KG*, Nydam, A*, Nott, Z., & Dux, PE 

### custom functions for running statistical analyses and plotting 
get.data <- function(data, subs, N){
  # given a dataframe, the list of subject numbers, & N (the number of subs required),
  # return a dataframe with that number of subs returned by replacement
  data[data$Subjects %in% sample(subs, size=N, replace=TRUE), ]
}

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