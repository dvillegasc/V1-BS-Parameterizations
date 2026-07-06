#' The Birnbaum-Saunders distribution - Second parameterization
#' 
#' @description
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the 
#' Birnbaum-Saunders distribution with 
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter (\code{mu > 0}).    
#' @param sigma parameter (\code{sigma > 0}).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).  
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{BS4}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{1}{2\sqrt{2\pi}} \left[ \frac{\sigma}{x\sqrt{x}} + \frac{\mu}{\sqrt{x}} \right] \exp\left( -\frac{1}{2} \left[ \frac{\sigma}{\sqrt{x}} - \mu\sqrt{x} \right]^2 \right)}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this 
#' parameterization 
#' \eqn{E(X) = \frac{\sigma \mu + 1/2}{\mu^2}} and 
#' \eqn{Var(X) = \frac{\sigma \mu + 5/4}{\mu^4}}.
#' 
#' @return 
#' \code{dBS4} gives the density, \code{pBS4} gives the distribution 
#' function, \code{qBS4} gives the quantile function, \code{rBS4}
#' generates random deviates and \code{hBS4} gives the hazard function.
#' 
#' @example examples/examples_dBS3.R
#' 
#' @export
dBS4 <- function(x, mu=1, sigma=0.5, log=FALSE){ #mu = mu   y  sigma = lambda
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS4 to BS (original)
  new_mu    <- (sigma/mu) #Beta
  new_sigma <- 1 / sqrt(mu*sigma) #Alfa
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS4
pBS4 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS4 to BS (original)
  new_mu    <- (sigma/mu)
  new_sigma <- 1 / sqrt(mu*sigma)
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS4
qBS4 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS4 to BS (original)
  new_mu    <- (sigma/mu)
  new_sigma <- 1 / sqrt(mu*sigma)
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS4
rBS4 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
      
      # Changing from BS4 to BS (original)
      new_mu    <- (sigma/mu)
      new_sigma <- 1 / sqrt(mu*sigma)
      
      r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
      r
}
#' @export
#' @rdname dBS4
hBS4 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
      
      h <- dBS4(x, mu, sigma) / pBS4(x, mu, sigma, lower.tail=FALSE)
      h
}

  