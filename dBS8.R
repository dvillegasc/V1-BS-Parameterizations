#' The Birnbaum-Saunders distribution - Sixth parameterization (Based on the variance 2)
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
#' @param mu parameter representing the variance \eqn{\sigma^2} (\code{mu > 0}).    
#' @param sigma parameter representing the shape \eqn{\gamma} (\code{sigma > 1}).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).  
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{BS8}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{\sqrt{\sigma}}{2\sqrt{2\pi\mu}} \left[ \left\{ \frac{1}{2x} \sqrt{\frac{5\mu}{\sigma(\sigma-1)}} \right\}^{1/2} + \left\{ \frac{1}{2x} \sqrt{\frac{5\mu}{\sigma(\sigma-1)}} \right\}^{3/2} \right] \exp\left( -\frac{5}{8(\sigma-1)} \left[ \frac{2x\sqrt{\sigma(\sigma-1)}}{\sqrt{5\mu}} + \frac{\sqrt{5\mu}}{2x\sqrt{\sigma(\sigma-1)}} - 2 \right] \right)}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>1}. In this parameterization, 
#' \eqn{E(X) = \frac{(2\sigma + 3)\sqrt{\mu}}{\sqrt{20\sigma(\sigma - 1)}}} and 
#' \eqn{Var(X) = \mu}.
#' 
#' @return 
#' \code{dBS8} gives the density, \code{pBS8} gives the distribution 
#' function, \code{qBS8} gives the quantile function, \code{rBS8}
#' generates random deviates and \code{hBS8} gives the hazard function.
#' 
#' @example examples/examples_dBS8.R
#' 
#' @export
dBS8 <- function(x, mu=1, sigma=0.5, log=FALSE){ #mu = varianza   y  sigma = alpha
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 1))  stop(paste("sigma must be positive", "\n", "")) #(based on the variance 2)
  
  # Changing from BS to BS8 (original)
  new_mu    <- (sqrt(5*mu))/(2 * sqrt(sigma * (sigma-1))) #Beta
  new_sigma <-  (2 * sqrt(sigma-1))/sqrt(5) #Alfa
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS8
pBS8 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 1))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS8 (original)
  new_mu    <- (sqrt(5*mu))/(2 * sqrt(sigma * (sigma-1)))
  new_sigma <-  (2 * sqrt(sigma-1))/sqrt(5)
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS8
qBS8 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 1)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS8 (original)
  new_mu    <- (sqrt(5*mu))/(2 * sqrt(sigma * (sigma-1)))
  new_sigma <-  (2 * sqrt(sigma-1))/sqrt(5)
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS8
rBS8 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 1))
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS to BS8 (original)
  new_mu    <- (sqrt(5*mu))/(2 * sqrt(sigma * (sigma-1)))
  new_sigma <-  (2 * sqrt(sigma-1))/sqrt(5)
  
  r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
  r
}
#' @export
#' @rdname dBS8
hBS8 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 1))
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dBS8(x, mu, sigma) / pBS8(x, mu, sigma, lower.tail=FALSE)
  h
}

