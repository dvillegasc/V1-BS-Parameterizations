#' The Birnbaum-Saunders family - Eleventh parameterization (Fourth Tweedie)
#' 
#' @description 
#' The function \code{BS13()} defines the Birnbaum-Saunders distribution, 
#' a two-parameter distribution, for a \code{gamlss.family} object 
#' to be used in GAMLSS fitting using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default 
#' for the mu parameter (representing the scale \eqn{\omega}).
#' @param sigma.link defines the sigma.link, with "log" link as the default 
#' for the sigma parameter (representing the shape \eqn{\psi}).
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Ahmed, S. E. (2012). 
#' On new parameterizations of the Birnbaum-Saunders distribution. 
#' Pakistan Journal of Statistics, 28(1), 1-26.
#' 
#' @seealso \link{dBS13}.
#' 
#' @details 
#' The Birnbaum-Saunders distribution with parameters \code{mu} and \code{sigma} 
#' (where \code{mu} represents \eqn{\omega} and \code{sigma} represents \eqn{\psi}) 
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{1}{\sqrt{2\pi}} \exp\left( -\frac{\sigma}{2} \left[ \frac{x\sigma}{\mu} + \frac{\mu}{x\sigma} - 2 \right] \right) \frac{[x\sigma + \mu]}{2\sqrt{\mu x^3}}}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization, 
#' \eqn{E(X) = \frac{\mu}{\sigma} + \frac{\mu}{2\sigma^2}} and 
#' \eqn{Var(X) = \frac{\mu^2}{\sigma^3} + \frac{5\mu^2}{4\sigma^4}}.
#' 
#' @returns Returns a \code{gamlss.family} object which can be used to fit a 
#' BS13 distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BS13.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS13 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS13", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "BS13", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("BS13", "Birnbaum-Saunders - Thirteenth parameterization"),
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
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           db_dm <- 1 / sigma
           da_dm <- 0 # <--
           
           term3 <- (1 / (y + b0)) * db_dm
           term4 <- (-1 / (2 * b0)) * db_dm
           term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           result <- term3 + term4 + term5
           return(result)
         },
         
         dldd = function(y, mu, sigma) { 
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           db_ds <- -b0 / sigma
           da_ds <- -1 / (2 * sigma * sqrt(sigma))
           
           term1 <- (-1 / a0) * da_ds
           term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           term3 <- (1 / (y + b0)) * db_ds
           term4 <- (-1 / (2 * b0)) * db_ds
           term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           
           result <- term1 + term2 + term3 + term4 + term5
           return(result)
         },
         
         # Second derivatives
         
         d2ldm2 = function(y, mu, sigma) {
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           db_dm <- 1 / sigma
           
           t3 <- (1 / (y + b0)) * db_dm
           t4 <- (-1 / (2 * b0)) * db_dm
           t5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           dldm <- t3 + t4 + t5
           
           return(-dldm * dldm) 
         },
         
         d2ldd2 = function(y, mu, sigma) {
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           
           db_ds <- -b0 / sigma
           da_ds <- -1 / (2 * sigma * sqrt(sigma))
           
           t1 <- (-1 / a0) * da_ds
           t2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           t3 <- (1 / (y + b0)) * db_ds
           t4 <- (-1 / (2 * b0)) * db_ds
           t5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           dldd <- t1 + t2 + t3 + t4 + t5
           
           return(-dldd * dldd)
         },
         
         d2ldmdd = function(y, mu, sigma) {
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           
           db_dm <- 1 / sigma
           da_ds <- -1 / (2 * sigma * sqrt(sigma))
           db_ds <- -b0 / sigma
           
           # dldm
           m3 <- (1 / (y + b0)) * db_dm
           m4 <- (-1 / (2 * b0)) * db_dm
           m5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           dldm <- m3 + m4 + m5
           
           # dldd
           d1 <- (-1 / a0) * da_ds
           d2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           d3 <- (1 / (y + b0)) * db_ds
           d4 <- (-1 / (2 * b0)) * db_ds
           d5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           dldd <- d1 + d2 + d3 + d4 + d5
           
           return(-dldm * dldd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS13(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS13", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(median(y), length(y))}),
         sigma.initial = expression({sigma <- rep(1, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
