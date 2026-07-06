#' The Birnbaum-Saunders distribution - Fourth parameterization (Based on the mean)
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
#' @param mu parameter representing the mean (\code{mu > 0}).    
#' @param sigma parameter representing the shape \eqn{\alpha} (\code{sigma > 0}).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).  
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{BS6}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{\exp(1/\sigma^2)\sqrt{2+\sigma^2}}{4\sigma\sqrt{\pi\mu}x^{3/2}} \left[ x + \frac{2\mu}{2+\sigma^2} \right] \exp\left( -\frac{1}{2\sigma^2} \left[ \frac{\{2+\sigma^2\}x}{2\mu} + \frac{2\mu}{\{2+\sigma^2\}x} \right] \right)}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization, 
#' \eqn{E(X) = \mu} and 
#' \eqn{Var(X) = [\mu\sigma]^2 \left[ \frac{4+5\sigma^2}{(2+\sigma^2)^2} \right]}.
#' 
#' @return 
#' \code{dBS6} gives the density, \code{pBS6} gives the distribution 
#' function, \code{qBS6} gives the quantile function, \code{rBS6}
#' generates random deviates and \code{hBS6} gives the hazard function.
#' 
#' @example examples/examples_dBS6.R
#' 
#' @export
dBS6 <- function(x, mu=1, sigma=0.5, log=FALSE){ #mu = mu   y  sigma = Alfa (original)
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS6 to BS (original)
  new_mu    <- (2 * mu)/(2 + sigma ^ 2) #Beta
  new_sigma <-  sigma #Alfa
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS6
pBS6 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS66 to BS6 (original)
  new_mu    <- (2 * mu)/(2 + sigma ^ 2)
  new_sigma <-  sigma
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS6
qBS6 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS66 to BS6 (original)
  new_mu    <- (2 * mu)/(2 + sigma ^ 2)
  new_sigma <-  sigma
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS6
rBS6 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS66 to BS6 (original)
  new_mu    <- (2 * mu)/(2 + sigma ^ 2)
  new_sigma <-  sigma
  
  r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
  r
}
#' @export
#' @rdname dBS6
hBS6 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dBS6(x, mu, sigma) / pBS6(x, mu, sigma, lower.tail=FALSE)
  h
}

