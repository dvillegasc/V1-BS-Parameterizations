
#----------------------- BS4 Standar ------------------------

dldm = function(y, mu, sigma) {
  a0 <- 1 / sqrt(mu * sigma)
  b0 <- sigma / mu
  da_dm <- -a0 / (2 * mu)
  db_dm <- -b0 / mu
  
  term1 <- (-1 / a0) * da_dm
  term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_dm
  term3 <- (1 / (y + b0)) * db_dm
  term4 <- (-1 / (2 * b0)) * db_dm
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  result <- term1 + term2 + term3 + term4 + term5
  return(result)
}

dldd = function(y, mu, sigma) { 
  a0 <- 1 / sqrt(mu * sigma)
  b0 <- sigma / mu
  da_ds <- -a0 / (2 * sigma)
  db_ds <- b0 / sigma

  term1 <- (-1 / a0) * da_ds
  term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  term3 <- (1 / (y + b0)) * db_ds
  term4 <- (-1 / (2 * b0)) * db_ds
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  
  result <- term1 + term2 + term3 + term4 + term5
  return(result)
}

# Second derivatives
d2ldm2 = function(y, mu, sigma) {
  a0 <- 1 / sqrt(mu * sigma)
  b0 <- sigma / mu
  da_dm <- -a0 / (2 * mu)
  db_dm <- -b0 / mu
  
  t1 <- (-1 / a0) * da_dm
  t2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_dm
  t3 <- (1 / (y + b0)) * db_dm
  t4 <- (-1 / (2 * b0)) * db_dm
  t5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  dldm <- t1 + t2 + t3 + t4 + t5
  
  return(-dldm * dldm) 
}

d2ldd2 = function(y, mu, sigma) {
  a0 <- 1 / sqrt(mu * sigma)
  b0 <- sigma / mu
  da_ds <- -a0 / (2 * sigma)
  db_ds <- b0 / sigma
  
  t1 <- (-1 / a0) * da_ds
  t2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  t3 <- (1 / (y + b0)) * db_ds
  t4 <- (-1 / (2 * b0)) * db_ds
  t5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  dldd <- t1 + t2 + t3 + t4 + t5
  
  return(-dldd * dldd)
}

d2ldmdd = function(y, mu, sigma) {
  a0 <- 1 / sqrt(mu * sigma)
  b0 <- sigma / mu
  
  da_dm <- -a0 / (2 * mu)
  db_dm <- -b0 / mu
  da_ds <- -a0 / (2 * sigma)
  db_ds <- b0 / sigma
  
  # dldm
  m1 <- (-1 / a0) * da_dm
  m2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_dm
  m3 <- (1 / (y + b0)) * db_dm
  m4 <- (-1 / (2 * b0)) * db_dm
  m5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  dldm <- m1 + m2 + m3 + m4 + m5
  
  # dldd
  d1 <- (-1 / a0) * da_ds
  d2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  d3 <- (1 / (y + b0)) * db_ds
  d4 <- (-1 / (2 * b0)) * db_ds
  d5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  dldd <- d1 + d2 + d3 + d4 + d5
  
  return(-dldm * dldd)
}

#----------------------- BS4 -------------------------------

# First derivatives
dldm = function(y, mu, sigma) {
  result <- (y / (sigma + mu * y)) + sigma - (mu * y)
  return(result)
}

dldd = function(y, sigma, mu) {
  result <- (1 / (sigma + mu * y)) + mu - (sigma / y)
  return(result)
}

# Second derivatives

d2ldm2 = function(y, sigma, mu) {
  result <- (y / (sigma + mu * y)) + sigma - (mu * y)
  return(-result * result)
}

d2ldd2 = function(y, sigma, mu) {
  result <- (1 / (sigma + mu * y)) + mu - (sigma / y)
  return(-result * result)
}

d2ldmdd = function(y, sigma, mu) {
  
  dldm <- (y / (sigma + mu * y)) + sigma - (mu * y)
  
  dldd <- (1 / (sigma + mu * y)) + mu - (sigma / y)
  
  d2ldmdd <- -dldm * dldd
  return(d2ldmdd)
}

#----------------------- BS5 -------------------------------

# First derivatives
dldm = function(y, mu, sigma) {
  
  term1 <- sigma / (y * (sigma + 1) + sigma * mu)
  term2 <- -1 / (2 * mu)
  term3 <- (y * (sigma + 1)) / (4 * mu^2)
  term4 <- -(sigma^2) / (4 * y * (sigma + 1))
  
  result <- (term1 + term2 + term3 + term4)
  return(result)
}

dldd = function(y, mu, sigma) { 
  
  term1 <- 0.5 + (1 / (2 * (sigma + 1)))
  term2 <- mu / ((y * (sigma + 1) + sigma * mu) * (sigma + 1))
  term3 <- -y / (4 * mu)
  term4 <- -(mu * (sigma^2 + 2 * sigma)) / (4 * y * (sigma + 1)^2)
  
  result <- (term1 + term2 + term3 + term4)
  return(result)
}

# Second derivatives

d2ldm2 = function(y, mu, sigma) {
  
  term1 <- sigma / (y * (sigma + 1) + sigma * mu)
  term2 <- -1 / (2 * mu)
  term3 <- (y * (sigma + 1)) / (4 * mu^2)
  term4 <- -(sigma^2) / (4 * y * (sigma + 1))
  dldm <- term1 + term2 + term3 + term4
  
  return(-dldm * dldm) 
}

d2ldd2 = function(y, mu, sigma) {
  term1 <- 0.5 + (1 / (2 * (sigma + 1)))
  term2 <- mu / ((y * (sigma + 1) + sigma * mu) * (sigma + 1))
  term3 <- -y / (4 * mu)
  term4 <- -(mu * (sigma^2 + 2 * sigma)) / (4 * y * (sigma + 1)^2)
  dldd <- term1 + term2 + term3 + term4
  
  return(-dldd * dldd)
}

d2ldmdd = function(y, mu, sigma) {
  
  # dldm
  dm1 <- sigma / (y * (sigma + 1) + sigma * mu)
  dm2 <- -1 / (2 * mu)
  dm3 <- (y * (sigma + 1)) / (4 * mu^2)
  dm4 <- -(sigma^2) / (4 * y * (sigma + 1))
  dldm <- dm1 + dm2 + dm3 + dm4
  
  # dldd
  dd1 <- 0.5 + (1 / (2 * (sigma + 1)))
  dd2 <- mu / ((y * (sigma + 1) + sigma * mu) * (sigma + 1))
  dd3 <- -y / (4 * mu)
  dd4 <- -(mu * (sigma^2 + 2 * sigma)) / (4 * y * (sigma + 1)^2)
  dldd <- dd1 + dd2 + dd3 + dd4
  
  return(-dldm * dldd)
}

#----------------------- BS6 -------------------------------

# First derivatives
dldm = function(y, mu, sigma) {
  b <- (2 * mu) / (2 + sigma^2)
  db_dm <- 2 / (2 + sigma^2)
  
  term1 <- (1 / (y + b)) * db_dm
  term2 <- -1 / (2 * mu)
  term3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  
  result <- term1 + term2 + term3
  return(result)
}

dldd = function(y, mu, sigma) { 
  b <- (2 * mu) / (2 + sigma^2)
  db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)
  
  term1 <- -1 / sigma
  term2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  term3 <- (1 / (y + b)) * db_ds
  term4 <- (-1 / (2 * b)) * db_ds
  term5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  
  result <- term1 + term2 + term3 + term4 + term5
  return(result)
}

# Second derivatives

d2ldm2 = function(y, mu, sigma) {
  b <- (2 * mu) / (2 + sigma^2)
  db_dm <- 2 / (2 + sigma^2)
  
  term1 <- (1 / (y + b)) * db_dm
  term2 <- -1 / (2 * mu)
  term3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  
  dldm <- term1 + term2 + term3
  
  return(-dldm * dldm) 
}

d2ldd2 = function(y, mu, sigma) {
  b <- (2 * mu) / (2 + sigma^2)
  db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)
  
  term1 <- -1 / sigma
  term2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  term3 <- (1 / (y + b)) * db_ds
  term4 <- (-1 / (2 * b)) * db_ds
  term5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  
  dldd <- term1 + term2 + term3 + term4 + term5
  
  return(-dldd * dldd)
}

d2ldmdd = function(y, mu, sigma) {
  b <- (2 * mu) / (2 + sigma^2)
  db_dm <- 2 / (2 + sigma^2)
  db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)
  
  # dldm
  m1 <- (1 / (y + b)) * db_dm
  m2 <- -1 / (2 * mu)
  m3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  dldm <- m1 + m2 + m3
  
  # dldd
  d1 <- -1 / sigma
  d2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  d3 <- (1 / (y + b)) * db_ds
  d4 <- (-1 / (2 * b)) * db_ds
  d5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  dldd <- d1 + d2 + d3 + d4 + d5
  
  return(-dldm * dldd)
}

#----------------------- BS7 -------------------------------

# First derivatives
dldm = function(y, mu, sigma) {
  b <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
  db_dm <- b / (2 * mu)
  
  term1 <- (1 / (y + b)) * db_dm
  term2 <- -1 / (2 * b) * db_dm
  term3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  
  result <- term1 + term2 + term3
  return(result)
}

dldd = function(y, mu, sigma) { 
  b <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
  num_deriv <- 4 + 10 * sigma^2
  den_deriv <- sigma * (4 + 5 * sigma^2)
  db_ds <- -b * (num_deriv / den_deriv)
  
  term1 <- -1 / sigma
  term2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  term3 <- (1 / (y + b)) * db_ds
  term4 <- -1 / (2 * b) * db_ds
  term5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  
  result <- term1 + term2 + term3 + term4 + term5
  return(result)
}

# Second derivatives

d2ldm2 = function(y, mu, sigma) {
  b <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
  db_dm <- b / (2 * mu)
  
  term1 <- (1 / (y + b)) * db_dm
  term2 <- -1 / (2 * b) * db_dm
  term3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  
  dldm <- term1 + term2 + term3
  
  return(-dldm * dldm) 
}

d2ldd2 = function(y, mu, sigma) {
  b <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
  
  num_deriv <- 4 + 10 * sigma^2
  den_deriv <- sigma * (4 + 5 * sigma^2)
  db_ds <- -b * (num_deriv / den_deriv)
  
  # dldd
  term1 <- -1 / sigma
  term2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  term3 <- (1 / (y + b)) * db_ds
  term4 <- -1 / (2 * b) * db_ds
  term5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  
  dldd <- term1 + term2 + term3 + term4 + term5
  
  return(-dldd * dldd)
}

d2ldmdd = function(y, mu, sigma) {
  b <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
  db_dm <- b / (2 * mu)
  num_deriv <- 4 + 10 * sigma^2
  den_deriv <- sigma * (4 + 5 * sigma^2)
  db_ds <- -b * (num_deriv / den_deriv)
  
  # dldm 
  m1 <- (1 / (y + b)) * db_dm
  m2 <- -1 / (2 * b) * db_dm
  m3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  dldm <- m1 + m2 + m3
  
  # dldd 
  d1 <- -1 / sigma
  d2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  d3 <- (1 / (y + b)) * db_ds
  d4 <- -1 / (2 * b) * db_ds
  d5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  dldd <- d1 + d2 + d3 + d4 + d5
  
  return(-dldm * dldd)
}





#----------------------- BS8 -------------------------------

# First derivatives
dldm = function(y, mu, sigma) {
  a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
  b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
  db_dm <- b0 / (2 * mu)
  
  term1 <- (1 / (y + b0)) * db_dm
  term2 <- -1 / (2 * b0) * db_dm
  term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  result <- term1 + term2 + term3
  return(result)
}

dldd = function(y, mu, sigma) { 
  a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
  b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
  
  da_ds <- 1 / sqrt(5 * (sigma - 1))
  db_ds <- -b0 * (2 * sigma - 1) / (2 * sigma * (sigma - 1))
  
  term1 <- (-1 / a0) * da_ds
  term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  term3 <- (1 / (y + b0)) * db_ds
  term4 <- (-1 / (2 * b0)) * db_ds
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  
  result <- term1 + term2 + term3 + term4 + term5
  return(result)
}

# Second derivatives

d2ldm2 = function(y, mu, sigma) {
  a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
  b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
  db_dm <- b0 / (2 * mu)
  
  term1 <- (1 / (y + b0)) * db_dm
  term2 <- -1 / (2 * b0) * db_dm
  term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  dldm <- term1 + term2 + term3
  
  return(-dldm * dldm) 
}

d2ldd2 = function(y, mu, sigma) {
  a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
  b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
  
  da_ds <- 1 / sqrt(5 * (sigma - 1))
  db_ds <- -b0 * (2 * sigma - 1) / (2 * sigma * (sigma - 1))
  
  term1 <- (-1 / a0) * da_ds
  term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  term3 <- (1 / (y + b0)) * db_ds
  term4 <- (-1 / (2 * b0)) * db_ds
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  
  dldd <- term1 + term2 + term3 + term4 + term5
  
  return(-dldd * dldd)
}

d2ldmdd = function(y, mu, sigma) {
  a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
  b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
  
  db_dm <- b0 / (2 * mu)
  da_ds <- 1 / sqrt(5 * (sigma - 1))
  db_ds <- -b0 * (2 * sigma - 1) / (2 * sigma * (sigma - 1))
  
  # dldm
  m1 <- (1 / (y + b0)) * db_dm
  m2 <- -1 / (2 * b0) * db_dm
  m3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  dldm <- m1 + m2 + m3
  
  # dldd
  d1 <- (-1 / a0) * da_ds
  d2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  d3 <- (1 / (y + b0)) * db_ds
  d4 <- (-1 / (2 * b0)) * db_ds
  d5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  dldd <- d1 + d2 + d3 + d4 + d5
  
  return(-dldm * dldd)
}

#----------------------- BS9 --------------------------------



#----------------------- BS10 -------------------------------



#----------------------- BS11 -------------------------------



#----------------------- BS12 -------------------------------



#----------------------- BS13 -------------------------------
