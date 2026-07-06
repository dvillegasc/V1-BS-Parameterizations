#' The Birnbaum-Saunders family - Santos-Neto et al. (2012) (P6 Based on the variance 2)
#' 
#' @author David Villegas Ceballos, \email{david.villegas1@udea.edu.co}
#' 
#' @description 
#' The function \code{BS8()} defines the Birnbaum-Saunders distribution, 
#' a two-parameter distribution, for a \code{gamlss.family} object 
#' to be used in GAMLSS fitting using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default 
#' for the mu parameter (representing the variance).
#' @param sigma.link defines the sigma.link, with "log" link as the default 
#' for the sigma parameter (representing the shape).
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{dBS8}.
#' 
#' @details 
#' The Birnbaum-Saunders distribution with parameters \code{mu} and \code{sigma} 
#' (where \code{mu} represents the true variance \eqn{\sigma^2} and \code{sigma} represents the shape parameter \eqn{\alpha}) 
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = 
#' \frac{\sqrt{\mu}}
#'      {2\sqrt{2\pi\sigma}}
#' \left[
#' \left\{
#' \frac{1}{2x}
#' \sqrt{\frac{5\sigma}{\mu(\mu-1)}}
#' \right\}^{1/2}
#' +
#' \left\{
#' \frac{1}{2x}
#' \sqrt{\frac{5\sigma}{\mu(\mu-1)}}
#' \right\}^{3/2}
#' \right]
#' \exp\left(
#' -\frac{5}{8(\mu-1)}
#' \left[
#' \frac{2x\sqrt{\mu(\mu-1)}}{\sqrt{5\sigma}}
#' +
#' \frac{\sqrt{5\sigma}}
#'      {2+\sqrt{\mu(\mu-1)}}
#' -2
#' \right]
#' \right)
#' }
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization, 
#' \eqn{E(X) = \frac{[2\mu+3]\sqrt{\sigma}}{\sqrt{20\mu(\mu-1)}}} and 
#' \eqn{Var(X) = \sigma}.
#' 
#' @returns Returns a \code{gamlss.family} object which can be used to fit a 
#' BS8 distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BS8.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS8 <- function(mu.link = "log", sigma.link = "log") {
  mstats <- checklink("mu.link", "BS8", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "BS8", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("BS8", "Birnbaum-Saunders - Seventh parameterization"),
         parameters = list(mu=TRUE, sigma=TRUE),
         nopar = 2,
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         
         # # First derivatives
         # dldm = function(y, mu, sigma) { 
         #   a0 <- (2 * sqrt(mu - 1)) / sqrt(5)
         #   b0 <- sqrt(5 * sigma) / (2 * sqrt(mu * (mu - 1)))
         #   
         #   da_ds <- 1 / sqrt(5 * (mu - 1))
         #   db_ds <- -b0 * (2 * mu - 1) / (2 * mu * (mu - 1))
         #   
         #   term1 <- (-1 / a0) * da_ds
         #   term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
         #   term3 <- (1 / (y + b0)) * db_ds
         #   term4 <- (-1 / (2 * b0)) * db_ds
         #   term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
         #   
         #   result <- term1 + term2 + term3 + term4 + term5
         #   return(result)
         # },
         # 
         # dldd = function(y, mu, sigma) {
         #   a0 <- (2 * sqrt(mu - 1)) / sqrt(5)
         #   b0 <- sqrt(5 * sigma) / (2 * sqrt(mu * (mu - 1)))
         #   db_dm <- b0 / (2 * sigma)
         #   
         #   term1 <- (1 / (y + b0)) * db_dm
         #   term2 <- -1 / (2 * b0) * db_dm
         #   term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
         #   
         #   result <- term1 + term2 + term3
         #   return(result)
         # },
         # 
         # # Second derivatives
         # d2ldm2 = function(y, mu, sigma) {
         #   a0 <- (2 * sqrt(mu - 1)) / sqrt(5)
         #   b0 <- sqrt(5 * sigma) / (2 * sqrt(mu * (mu - 1)))
         #   
         #   da_ds <- 1 / sqrt(5 * (mu - 1))
         #   db_ds <- -b0 * (2 * mu - 1) / (2 * mu * (mu - 1))
         #   
         #   term1 <- (-1 / a0) * da_ds
         #   term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
         #   term3 <- (1 / (y + b0)) * db_ds
         #   term4 <- (-1 / (2 * b0)) * db_ds
         #   term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
         #   
         #   dldd <- term1 + term2 + term3 + term4 + term5
         #   
         #   return(-dldd * dldd)
         # },
         # 
         # d2ldd2 = function(y, mu, sigma) {
         #   a0 <- (2 * sqrt(mu - 1)) / sqrt(5)
         #   b0 <- sqrt(5 * sigma) / (2 * sqrt(mu * (mu - 1)))
         #   db_dm <- b0 / (2 * sigma)
         #   
         #   term1 <- (1 / (y + b0)) * db_dm
         #   term2 <- -1 / (2 * b0) * db_dm
         #   term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
         #   
         #   dldm <- term1 + term2 + term3
         #   
         #   return(-dldm * dldm) 
         # },
         # 
         # d2ldmdd = function(y, mu, sigma) {
         #   a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
         #   b0 <- sqrt(5 * sigma) / (2 * sqrt(mu * (mu - 1)))
         #   
         #   db_dm <- b0 / (2 * sigma)
         #   da_ds <- 1 / sqrt(5 * (mu - 1))
         #   db_ds <- -b0 * (2 * mu - 1) / (2 * mu * (mu - 1))
         #   
         #   # dldm
         #   m1 <- (1 / (y + b0)) * db_dm
         #   m2 <- -1 / (2 * b0) * db_dm
         #   m3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
         #   dldm <- m1 + m2 + m3
         #   
         #   # dldd
         #   d1 <- (-1 / a0) * da_ds
         #   d2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
         #   d3 <- (1 / (y + b0)) * db_ds
         #   d4 <- (-1 / (2 * b0)) * db_ds
         #   d5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
         #   dldd <- d1 + d2 + d3 + d4 + d5
         #   
         #   return(-dldm * dldd)
         # },
         
         # First derivates
         
         dldm = function(y, mu, sigma) {
           dm   <- gamlss::numeric.deriv(dBS8(y, mu, sigma, log=TRUE),
                                         theta="mu",
                                         delta=0.00001)
           dldm <- as.vector(attr(dm, "gradient"))
           dldm
         },
         
         dldd = function(y, mu, sigma) {
           dd   <- gamlss::numeric.deriv(dBS8(y, mu, sigma, log=TRUE),
                                         theta="sigma",
                                         delta=0.00001)
           dldd <- as.vector(attr(dd, "gradient"))
           dldd
         },
         
         # Second derivates
         
         d2ldm2 = function(y, mu, sigma) {
           dm   <- gamlss::numeric.deriv(dBS8(y, mu, sigma, log=TRUE),
                                         theta="mu",
                                         delta=0.00001)
           dldm <- as.vector(attr(dm, "gradient"))
           d2ldm2 <- - dldm * dldm
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
           d2ldm2
         },
         
         d2ldmdd = function(y, mu, sigma) {
           dm   <- gamlss::numeric.deriv(dBS8(y, mu, sigma, log=TRUE),
                                         theta="mu",
                                         delta=0.00001)
           dldm <- as.vector(attr(dm, "gradient"))
           dd   <- gamlss::numeric.deriv(dBS8(y, mu, sigma, log=TRUE),
                                         theta="sigma",
                                         delta=0.00001)
           dldd <- as.vector(attr(dd, "gradient"))
           
           d2ldmdd <- - dldm * dldd
           d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
           d2ldmdd
         },
         
         d2ldd2  = function(y, mu, sigma) {
           dd   <- gamlss::numeric.deriv(dBS8(y, mu, sigma, log=TRUE),
                                         theta="sigma",
                                         delta=0.00001)
           dldd <- as.vector(attr(dd, "gradient"))
           d2ldd2 <- - dldd * dldd
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
           d2ldd2
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS8(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS8", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial = expression({mu <- rep((5*mean(y)^2+3*var(y)+5*mean(y)*sqrt(mean(y)^2+3*var(y)))/(2*(5*mean(y)^2-var(y))), length(y)) }),
         sigma.initial = expression({sigma <- rep(var(y), length(y))}),
         
         mu.valid = function(mu) all(mu > 1) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

