#' The Birnbaum-Saunders family - Second parameterization (Ahmed et al., 2008)
#'
#' @author David Villegas Ceballos, \email{david.villegas1@udea.edu.co}
#'
#' @description
#' The function \code{BS4()} defines the Birnbaum-Saunders distribution,
#' a two-parameter distribution, for a \code{gamlss.family} object
#' to be used in GAMLSS fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default
#' for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default
#' for the sigma parameter.
#'
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012).
#' On new parameterizations of the Birnbaum-Saunders distribution.
#' Pakistan Journal of Statistics, 28(1), 1-26.
#'
#' Ahmed, S. E., et al. (2008). Inference in an applied accelerated life test model
#' based on the Birnbaum-Saunders distribution. Journal of Statistical Computation
#' and Simulation, 78(9), 809-823.
#'
#' @seealso \link{dBS4}.
#'
#' @details
#' The Birnbaum-Saunders distribution with parameters \code{mu} and \code{sigma}
#' has density given by
#'
#' \eqn{f(x|\mu,\sigma) = \frac{1}{2\sqrt{2\pi}} \left[ \frac{\sigma}{x\sqrt{x}} + \frac{\mu}{\sqrt{x}} \right] \exp\left( -\frac{1}{2} \left[ \frac{\sigma}{\sqrt{x}} - \mu\sqrt{x} \right]^2 \right)}
#'
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization,
#' \eqn{E(X) = \frac{\sigma \mu + 1/2}{\mu^2}} and
#' \eqn{Var(X) = \frac{\sigma \mu + 5/4}{\mu^4}}.
#'
#' @returns Returns a \code{gamlss.family} object which can be used to fit a
#' BS4 distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_BS4.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS4 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS4", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "BS4", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("BS4", "Birnbaum-Saunders - Fourth parameterization"),
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
         
         # First derivatives
         dldm = function(y, mu, sigma) {
           result <- (y / (sigma + mu * y)) + sigma - (mu * y)
           return(result)
         },
         
         dldd = function(y, sigma, mu) {
           result <- (1 / (sigma + mu * y)) + mu - (sigma / y)
           return(result)
         },
         
         # Second derivatives
         
         d2ldm2 = function(y, sigma, mu) {
           result <- (y / (sigma + mu * y)) + sigma - (mu * y)
           return(-result * result)
         },
         
         d2ldd2 = function(y, sigma, mu) {
           result <- (1 / (sigma + mu * y)) + mu - (sigma / y)
           return(-result * result)
         },
         
         d2ldmdd = function(y, sigma, mu) {
           
           dldm <- (y / (sigma + mu * y)) + sigma - (mu * y)
           
           dldd <- (1 / (sigma + mu * y)) + mu - (sigma / y)
           
           d2ldmdd <- -dldm * dldd
           return(d2ldmdd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS4(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS4", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep((1 / sqrt(median(y))), length(y))}),
         sigma.initial = expression({sigma <- rep(sqrt(median(y)), length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
