#' The Birnbaum-Saunders distribution - Eleventh parameterization (Fourth Tweedie)
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
#' @param mu parameter representing \eqn{\omega} (\code{mu > 0}).    
#' @param sigma parameter representing \eqn{\psi} (\code{sigma > 0}).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).  
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{BS13}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{1}{\sqrt{2\pi}} \exp\left( -\frac{\sigma}{2} \left[ \frac{x\sigma}{\mu} + \frac{\mu}{x\sigma} - 2 \right] \right) \frac{[x\sigma + \mu]}{2\sqrt{\mu x^3}}}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization, 
#' \eqn{E(X) = \frac{\mu}{\sigma} + \frac{\mu}{2\sigma^2}} and 
#' \eqn{Var(X) = \frac{\mu^2}{\sigma^3} + \frac{5\mu^2}{4\sigma^4}}.
#' 
#' @return 
#' \code{dBS13} gives the density, \code{pBS13} gives the distribution 
#' function, \code{qBS13} gives the quantile function, \code{rBS13}
#' generates random deviates and \code{hBS13} gives the hazard function.
#' 
#' @example examples/examples_dBS13.R
#' 
#' @export
dBS13 <- function(x, mu=1, sigma=0.5, log=FALSE){ #mu = omega   y  sigma = Ψ
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) #(based on the variance 2)
  
  # Changing from BS to BS13 (original)
  new_mu    <- mu / sigma #Beta
  new_sigma <-  1 / sqrt(sigma) #Alfa
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS13
pBS13 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS13 (original)
  new_mu    <- mu / sigma
  new_sigma <-  1 / sqrt(sigma)
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS13
qBS13 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS13 (original)
  new_mu    <- mu / sigma
  new_sigma <-  1 / sqrt(sigma)  
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS13
rBS13 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS13 (original)
  new_mu    <- mu / sigma
  new_sigma <-  1 / sqrt(sigma)
  
  r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
  r
}
#' @export
#' @rdname dBS13
hBS13 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dBS13(x, mu, sigma) / pBS13(x, mu, sigma, lower.tail=FALSE)
  h
}

