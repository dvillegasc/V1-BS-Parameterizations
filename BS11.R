#' The Birnbaum-Saunders family - Ninth parameterization (Second Tweedie)
#' 
#' @description 
#' The function \code{BS11()} defines the Birnbaum-Saunders distribution, 
#' a two-parameter distribution, for a \code{gamlss.family} object 
#' to be used in GAMLSS fitting using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default 
#' for the mu parameter (representing \eqn{\beta}).
#' @param sigma.link defines the sigma.link, with "log" link as the default 
#' for the sigma parameter (representing \eqn{\omega}).
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{dBS11}.
#' 
#' @details 
#' The Birnbaum-Saunders distribution with parameters \code{mu} and \code{sigma} 
#' (where \code{mu} represents \eqn{\beta} and \code{sigma} represents \eqn{\omega}) 
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{1}{\sqrt{2\pi}} \exp\left( -\frac{\sigma}{2\mu} \left[ \frac{x}{\mu} + \frac{\mu}{x} - 2 \right] \right) \frac{[x + \mu]\sqrt{\sigma}}{2\mu\sqrt{x^3}}}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization, 
#' \eqn{E(X) = \mu + \frac{\mu^2}{2\sigma}} and 
#' \eqn{Var(X) = \frac{\mu^3}{\sigma} + \frac{5\mu^4}{4\sigma^2}}.
#' 
#' @returns Returns a \code{gamlss.family} object which can be used to fit a 
#' BS11 distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BS11.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS11 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS11", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "BS11", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("BS11", "Birnbaum-Saunders - Eleventh parameterization"),
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
           b0 <- mu
           a0 <- sqrt(mu / sigma)
           db_dm <- 1
           da_dm <- 1 / (2 * sqrt(mu * sigma))
           
           term1 <- (-1 / a0) * da_dm
           term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_dm
           
           term3 <- (1 / (y + b0)) * db_dm
           term4 <- (-1 / (2 * b0)) * db_dm
           term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           result <- term1 + term2 + term3 + term4 + term5
           return(result)
         },
         
         dldd = function(y, mu, sigma) { 
           b0 <- mu
           a0 <- sqrt(mu / sigma)
           db_ds <- 0 # <--
           da_ds <- -a0 / (2 * sigma)
           
           term1 <- (-1 / a0) * da_ds
           term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           
           result <- term1 + term2
           return(result)
         },
         
         # Second derivatives
         
         d2ldm2 = function(y, mu, sigma) {
           b0 <- mu
           a0 <- sqrt(mu / sigma)
           
           db_dm <- 1
           da_dm <- 1 / (2 * sqrt(mu * sigma))
           
           # dldm
           t1 <- (-1 / a0) * da_dm
           t2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_dm
           t3 <- (1 / (y + b0)) * db_dm
           t4 <- (-1 / (2 * b0)) * db_dm
           t5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           dldm <- t1 + t2 + t3 + t4 + t5
           
           return(-dldm * dldm) 
         },
         
         d2ldd2 = function(y, mu, sigma) {
           b0 <- mu
           a0 <- sqrt(mu / sigma)
           
           da_ds <- -a0 / (2 * sigma)
           db_ds <- 0 # <--
           
           t1 <- (-1 / a0) * da_ds
           t2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           dldd <- t1 + t2
           
           return(-dldd * dldd)
         },
         
         d2ldmdd = function(y, mu, sigma) {
           b0 <- mu
           a0 <- sqrt(mu / sigma)
           
           db_dm <- 1
           da_dm <- 1 / (2 * sqrt(mu * sigma))
           da_ds <- -a0 / (2 * sigma)
           db_ds <- 0 # <--
           
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
           dldd <- d1 + d2
           
           return(-dldm * dldd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS11(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS11", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(median(y), length(y))}),
         sigma.initial = expression({sigma <- rep(1, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
