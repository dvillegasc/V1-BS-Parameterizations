#' The Birnbaum-Saunders distribution - Ninth parameterization (Second Tweedie)
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
#' @param mu parameter representing \eqn{\beta} (\code{mu > 0}).    
#' @param sigma parameter representing \eqn{\omega} (\code{sigma > 0}).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).  
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{BS11}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{1}{\sqrt{2\pi}} \exp\left( -\frac{\sigma}{2\mu} \left[ \frac{x}{\mu} + \frac{\mu}{x} - 2 \right] \right) \frac{[x + \mu]\sqrt{\sigma}}{2\mu\sqrt{x^3}}}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization, 
#' \eqn{E(X) = \mu + \frac{\mu^2}{2\sigma}} and 
#' \eqn{Var(X) = \frac{\mu^3}{\sigma} + \frac{5\mu^4}{4\sigma^2}}.
#' 
#' @return 
#' \code{dBS11} gives the density, \code{pBS11} gives the distribution 
#' function, \code{qBS11} gives the quantile function, \code{rBS11}
#' generates random deviates and \code{hBS11} gives the hazard function.
#' 
#' @example examples/examples_dBS11.R
#' 
#' @export
dBS11 <- function(x, mu=1, sigma=0.5, log=FALSE){ #mu = β   y  sigma = ώ
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) #(based on the variance 2)
  
  # Changing from BS to BS11 (original)
  new_mu    <- mu #Beta
  new_sigma <-  sqrt(mu / sigma) #Alfa
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS11
pBS11 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS11 (original)
  new_mu    <- mu 
  new_sigma <-  sqrt(mu / sigma) 
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS11
qBS11 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS11 (original)
  new_mu    <- mu 
  new_sigma <-  sqrt(mu / sigma)  
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS11
rBS11 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS11 (original)
  new_mu    <- mu 
  new_sigma <-  sqrt(mu / sigma) 
  
  r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
  r
}
#' @export
#' @rdname dBS11
hBS11 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dBS11(x, mu, sigma) / pBS11(x, mu, sigma, lower.tail=FALSE)
  h
}

