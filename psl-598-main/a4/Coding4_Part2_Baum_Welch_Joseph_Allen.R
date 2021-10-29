library(purrr)

forward_prob = function(x, params) {
  # OUtput the forward probability matrix alp
  ## alp: t-by-mz, (t, i) entry = P(x_{1 : t}, Z_t = i)
  t   = length(x)
  mz  = params$mz
  a   = params$A
  b   = params$B
  w   = params$w
  alp = matrix(0, t, mz)
  
  # fill in the first row of alp
  alp[1, ] = w * b[, x[1]]
  
  # recursively compute the remaining rows of alp
  for (t in 2 : t) {
    temp = alp[t - 1, ] %*% a
    alp[t, ] = temp * b[, x[t]]
  }
  return(alp)
}

backward_prob = function(x, params) {
  # Output the backward probability matrix beta
  ## beta: t-by-mz, (t, i) entry = P(x_{1 : t}, Z_t = i)
  t    = length(x)
  mz   = params$mz
  a    = params$A
  b    = params$B
  w    = params$w
  beta = matrix(1, t, mz)
  
  # the last row of beta is all 1.
  # recursively compute the previous rows of beta
  for (t in (t - 1) : 1) {
    temp = as.matrix(beta[t + 1, ] * b[, x[t + 1]])   # make temp a column vector
    beta[t, ] = t(a %*% temp)
  }
  return(beta)
}

bw_onestep = function(x, params) {
  # Input:
  # x: t-by-1 observation sequence
  # para: mx, mz, and current paarameter values for 
  ##  A: initial estimate for mz-by-mz transition matrix
  ##  B: initial estimate for mz-by-mx emission matrix
  ##  w: initial estimate for mz-by-1 initial distribution over Z_1
  # Output the updated parameters after one iteration
  # We DO NOT update the initial distribution w
  
  t    = length(x)
  mz   = params$mz
  mx   = params$mx
  a    = params$A
  b    = params$B
  w    = params$w
  alp  = forward_prob(x, params)
  beta = backward_prob(x, params)
  
  my_gamma = array(0, dim = c(mz, mz, t - 1))
  
  # Compute Gamma
  
  for (t in 1 : t - 1) {
    for (i in 1 : mz) {
      for (j in 1 : mz) {
        my_gamma[i, j, t] = alp[t, i] * beta[t + 1, j] * b[j, x[t + 1]] * a[i, j]
      }
    }
  }
  
  # M-step for parameter A
  
  # a = rowSums(my_gamma, dims = 2)
  a = rowSums(my_gamma, dims = 2) / rowSums(rowSums(my_gamma, dims = 2))
  
  # M-step for parameter B
  
  temp = cbind(apply(my_gamma, c(1, 3), sum), colSums(my_gamma[, , t - 1]))
  
  for (i in 1 : mx) {
    b[, i] = rowSums(temp[, which(x == i)])
  }
  
  b = b / rowSums(b)
  
  params$A = a
  params$B = b
  
  return(params)
}

my_bw = function(x, params) {
  # INPUT:
  ## x: t-by-1 observation sequence
  ## params: initial parameter value
  ## Output updated params value (A and B; we do not update w)
  
  for (i in 1 : 100) {
    params = bw_onestep(x, params)
  }
  
  return (params)
  
}

data = scan("coding4_part2.txt")

mz = 2
mx = 3
ini.w = rep(1, mz); ini.w = ini.w / sum(ini.w)
ini.A = matrix(1, 2, 2); ini.A = ini.A / rowSums(ini.A)
ini.B = matrix(1:6, 2, 3); ini.B = ini.B / rowSums(ini.B)
ini.para = list(mz = 2, mx = 3, w = ini.w,
                A = ini.A, B = ini.B)

myout = my_bw(data, ini.para)

