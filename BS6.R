#' The Birnbaum-Saunders family - Fourth parameterization (Based on the mean)
#' 
#' @description 
#' The function \code{BS6()} defines the Birnbaum-Saunders distribution, 
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
#' @seealso \link{dBS6}.
#' 
#' @details 
#' The Birnbaum-Saunders distribution with parameters \code{mu} and \code{sigma} 
#' (where \code{mu} represents the true mean and \code{sigma} represents the shape parameter \eqn{\alpha}) 
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{\exp(1/\sigma^2)\sqrt{2+\sigma^2}}{4\sigma\sqrt{\pi\mu}x^{3/2}} \left[ x + \frac{2\mu}{2+\sigma^2} \right] \exp\left( -\frac{1}{2\sigma^2} \left[ \frac{\{2+\sigma^2\}x}{2\mu} + \frac{2\mu}{\{2+\sigma^2\}x} \right] \right)}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this parameterization, 
#' \eqn{E(X) = \mu} and 
#' \eqn{Var(X) = [\mu\sigma]^2 \left[ \frac{4+5\sigma^2}{(2+\sigma^2)^2} \right]}.
#' 
#' @returns Returns a \code{gamlss.family} object which can be used to fit a 
#' BS6 distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BS6.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS6 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS6", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "BS6", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("BS6", "Birnbaum-Saunders - Sixth parameterization"),
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
           a0 <- sigma
           b0 <- (2 * mu) / (2 + sigma^2)
           db_dm <- 2 / (2 + sigma^2)
           
           term1 <- (1 / (y + b0)) * db_dm
           term2 <- -1 / (2 * b0) * db_dm
           term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           result <- term1 + term2 + term3
           return(result)
         },
         
         dldd = function(y, mu, sigma) { 
           a0 <- sigma
           b0 <- (2 * mu) / (2 + sigma^2)
           da_ds <- 1
           db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)

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
           b0 <- (2 * mu) / (2 + sigma^2)
           db_dm <- 2 / (2 + sigma^2)
           
           term1 <- (1 / (y + b0)) * db_dm
           term2 <- -1 / (2 * b0) * db_dm
           term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           dldm <- term1 + term2 + term3
           return(-dldm * dldm) 
         },
         
         d2ldd2 = function(y, mu, sigma) {
           a0 <- sigma
           b0 <- (2 * mu) / (2 + sigma^2)
           
           da_ds <- 1
           db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)
           
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
           b0 <- (2 * mu) / (2 + sigma^2)
           
           db_dm <- 2 / (2 + sigma^2)
           da_ds <- 1
           db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)
           
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
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS6(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS6", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(mean(y), length(y))}),
         sigma.initial = expression({sigma <- rep(0.5, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
