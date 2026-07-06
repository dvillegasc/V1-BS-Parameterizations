#' The Birnbaum-Saunders distribution - Third parameterization (Based on GLM)
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
#' @param sigma precision parameter \eqn{\delta} (\code{sigma > 0}).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).  
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{BS5}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{\exp(\sigma/2)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}x^{3/2}} \left[ x + \frac{\sigma\mu}{\sigma+1} \right] \exp\left( -\frac{\sigma}{4} \left[ \frac{x(\sigma+1)}{\sigma\mu} + \frac{\sigma\mu}{x(\sigma+1)} \right] \right)}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this 
#' parameterization 
#' \eqn{E(X) = \mu} and 
#' \eqn{Var(X) = \mu^2 \left[ \frac{2\sigma+5}{(\sigma+1)^2} \right]}.
#' 
#' @return 
#' \code{dBS5} gives the density, \code{pBS5} gives the distribution 
#' function, \code{qBS5} gives the quantile function, \code{rBS5}
#' generates random deviates and \code{hBS5} gives the hazard function.
#' 
#' @example examples/examples_dBS5.R
#' 
#' @export
dBS5 <- function(x, mu=1, sigma=0.5, log=FALSE){ #mu = mu   y  sigma = precision
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS5 to BS (original)
  new_mu    <- (sigma * mu)/(sigma + 1) #Beta
  new_sigma <-  sqrt(2 / sigma) #Alfa
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS5
pBS5 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS5 to BS (original)
  new_mu    <- (sigma * mu)/(sigma + 1)
  new_sigma <-  sqrt(2 / sigma) 
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS5
qBS5 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS5 to BS (original)
  new_mu    <- (sigma * mu)/(sigma + 1)
  new_sigma <-  sqrt(2 / sigma) 
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS5
rBS5 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS5 to BS (original)
  new_mu    <- (sigma * mu)/(sigma + 1)
  new_sigma <-  sqrt(2 / sigma) 
  
  r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
  r
}
#' @export
#' @rdname dBS5
hBS5 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dBS5(x, mu, sigma) / pBS5(x, mu, sigma, lower.tail=FALSE)
  h
}

