#' The Birnbaum-Saunders distribution - Seventh parameterization (Bounded Var)
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
#' @param sigma parameter representing \eqn{\phi} (\code{sigma > 1}).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).  
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{BS9}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{\exp(1/[2\{\sigma-1\}])[x\sigma + \mu]}{4\sqrt{\pi\sigma[\sigma-1]\mu x^{3/2}}} \exp\left( -\frac{1}{4[\sigma-1]} \left[ \frac{x\sigma}{\mu} + \frac{\mu}{x\sigma} \right] \right)}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>1}. In this parameterization, 
#' \eqn{E(X) = \mu} and 
#' \eqn{Var(X) = \mu^2 \frac{[\sigma-1][5\sigma-3]}{\sigma^2}}.
#' 
#' @return 
#' \code{dBS9} gives the density, \code{pBS9} gives the distribution 
#' function, \code{qBS9} gives the quantile function, \code{rBS9}
#' generates random deviates and \code{hBS9} gives the hazard function.
#' 
#' @example examples/examples_dBS9.R
#' 
#' @export
dBS9 <- function(x, mu=1, sigma=0.5, log=FALSE){ #mu = μ   y  sigma = Φ
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 1))  stop(paste("sigma must be positive", "\n", "")) #(based on the variance 2)
  
  # Changing from BS to BS9 (original)
  new_mu    <- mu / sigma #Beta
  new_sigma <-  sqrt(2 * (sigma-1)) #Alfa
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS9
pBS9 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 1))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS9 (original)
  new_mu    <- mu / sigma 
  new_sigma <-  sqrt(2 * (sigma-1)) 
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS9
qBS9 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 1)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS9 (original)
  new_mu    <- mu / sigma 
  new_sigma <-  sqrt(2 * (sigma-1))
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS9
rBS9 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 1))
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS9 (original)
  new_mu    <- mu / sigma 
  new_sigma <-  sqrt(2 * (sigma-1))
  
  r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
  r
}
#' @export
#' @rdname dBS9
hBS9 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 1))
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dBS9(x, mu, sigma) / pBS9(x, mu, sigma, lower.tail=FALSE)
  h
}

