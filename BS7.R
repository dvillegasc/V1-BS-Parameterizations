#' The Birnbaum-Saunders family - Fifth parameterization (Based on the variance)
#' 
#' @description 
#' The function \code{BS7()} defines the Birnbaum-Saunders distribution, 
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
#' @seealso \link{dBS7}.
#' 
#' @details 
#' The Birnbaum-Saunders distribution with parameters \code{mu} and \code{sigma} 
#' (where \code{mu} represents the true variance \eqn{\sigma^2} and \code{sigma} represents the shape parameter \eqn{\alpha}) 
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{1}{\sqrt{2\pi}} \exp\left( -\frac{1}{2\sigma^2} \left[ \frac{\sigma\sqrt{4+5\sigma^2}}{2\sqrt{\mu}x^{-1}} + \frac{2\sqrt{\mu}\{x\sigma\}^{-1}}{\sqrt{4+5\sigma^2}} - 2 \right] \right) \times \left[ \frac{\{x\sigma\}^{-1/2}\{4+5\sigma^2\}^{1/4}}{2^{3/2}\mu^{1/4}} + \frac{\mu^{1/4}}{\{x\sigma\}^{3/2}\sqrt{2}\{4+5\sigma^2\}^{1/4}} \right]}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization, 
#' \eqn{E(X) = \frac{[2+\sigma^2]\sqrt{\mu}}{\sigma\sqrt{4+5\sigma^2}}} and 
#' \eqn{Var(X) = \mu}.
#' 
#' @returns Returns a \code{gamlss.family} object which can be used to fit a 
#' BS7 distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BS7.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS7 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS7", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "BS7", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("BS7", "Birnbaum-Saunders - Seventh parameterization"),
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
         
         # First derivates
         dldm = function(y, mu, sigma) {
           a0 <- sigma
           b0 <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
           db_dm <- b0 / (2 * mu)
           
           term1 <- (1 / (y + b0)) * db_dm
           term2 <- -1 / (2 * b0) * db_dm
           term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           result <- term1 + term2 + term3
           return(result)
         },
         
         dldd = function(y, mu, sigma) { 
           a0 <- sigma
           b0 <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
           da_ds <- 1
           db_ds <- -b0 * ((4 + 10 * sigma^2) / (sigma * (4 + 5 * sigma^2)))
           
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
           a0 <- sigma
           b0 <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
           db_dm <- b0 / (2 * mu)
           
           term1 <- (1 / (y + b0)) * db_dm
           term2 <- -1 / (2 * b0) * db_dm
           term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           dldm <- term1 + term2 + term3
           return(-dldm * dldm) 
         },
         
         d2ldd2 = function(y, mu, sigma) {
           a0 <- sigma
           b0 <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
           
           da_ds <- 1
           db_ds <- -b0 * ((4 + 10 * sigma^2) / (sigma * (4 + 5 * sigma^2)))
           
           term1 <- (-1 / a0) * da_ds
           term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           term3 <- (1 / (y + b0)) * db_ds
           term4 <- (-1 / (2 * b0)) * db_ds
           term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           
           dldd <- term1 + term2 + term3 + term4 + term5
           return(-dldd * dldd)
         },
         
         d2ldmdd = function(y, mu, sigma) {
           a0 <- sigma
           b0 <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
           
           db_dm <- b0 / (2 * mu)
           da_ds <- 1
           db_ds <- -b0 * ((4 + 10 * sigma^2) / (sigma * (4 + 5 * sigma^2)))
           
           # dldm 
           m1 <- (1 / (y + b0)) * db_dm
           m2 <- -1 / (2 * b0) * db_dm
           m3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           dldm <- m1 + m2 + m3
           
           # dldd 
           d1 <- (-1 / a0) * da_ds
           d2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           d3 <- (1 / (y + b0)) * db_ds
           d4 <- (-1 / (2 * b0)) * db_ds
           d5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           dldd <- d1 + d2 + d3 + d4 + d5
           
           return(-dldm * dldd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS7(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS7", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(var(y), length(y))}),
         sigma.initial = expression({sigma <- rep(0.5, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
