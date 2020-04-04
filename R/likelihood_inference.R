#' Approximate profile likelihood
#'
get.likelihood <- function(params, data, params_init = NULL, L = 53.5) {

    B <- data$B
    E <- data$E
    S <- data$S

    param.names <- c("pi", "rho", "r", "ip_mean", "ip_sd")

    if (all(param.names %in% names(params))) { # likelihood

        pi <- params["pi"]
        rho <- params["rho"]
        r <- params["r"]
        alpha <- (params["ip_mean"] / params["ip_sd"])^2
        beta <- params["ip_mean"] / params["ip_sd"]^2

        if (1 - pi + pi * rho * (1 - 2 / r / L) < 0) {
            return(-Inf)
        }

        n <- length(B)

        ll <- 2 * n * log(r) + n * alpha * log(beta / (beta + r)) - n * log((1 - pi + pi * rho * (1 - 2 / r / L)))
        for (i in 1:n) {
            if (B[i] == 0) {
                ll <- ll + log(1 - pi)
            } else {
                ll <- ll + log(pi / L * rho)
            }
            ll <- ll + r * (S[i] - L)
            ll <- ll + log(pgamma(S[i] - B[i], alpha, beta + r) - pgamma(max(0, S[i] - E[i]), alpha, beta + r))
        }

        return(ll)
    } else { # profile likelihood

        params_init <- params_init[setdiff(param.names, names(params))]

        if (!all(param.names %in% names(c(params, params_init)))) {
            stop("To compute profile likelihood, params and params_init should together contain the following entries: pi, rho, r, alpha, beta.")
        }

        optim(params_init,
              function(params_now) get.likelihood(c(params, params_now), data),
              control=list(fnscale=-1))$value

    }

}

#' Likelihood inference
#'
#'
#' @import rootSolve
likelihood.inference <- function(data) {

    params_init <- c(pi = 0.5,
                     rho = 2,
                     r = 0.2,
                     ip_mean = 6,
                     ip_sd = 2)

    suppressWarnings(tmp <- optim(params_init,
                                  function(params) get.likelihood(params, data),
                                  control=list(fnscale=-1)))
    params_mle <- tmp$par
    ml <- tmp$value

    pi_mle <- params_mle["pi"]
    rho_mle <- params_mle["rho"]

    r_mle <- params_mle["r"]
    r_lrt <- function(r_now) ml - get.likelihood(c(r = r_now), data, params_mle) - qchisq(0.95, 1) / 2
    suppressWarnings(r_lower <- tryCatch(uniroot(r_lrt, c(0.1, r_mle))$root, error = function(e) NA))
    suppressWarnings(r_upper <- tryCatch(uniroot(r_lrt, c(r_mle, 0.5))$root, error = function(e) NA))

    ip_mean_mle <- params_mle["ip_mean"]
    ip_mean_lrt <- function(ip_mean_now) ml - get.likelihood(c(ip_mean = ip_mean_now), data, params_mle) - qchisq(0.95, 1) / 2
    suppressWarnings(ip_mean_lower <- tryCatch(uniroot(ip_mean_lrt, c(2, ip_mean_mle))$root, error = function(e) NA))
    suppressWarnings(ip_mean_upper <- tryCatch(uniroot(ip_mean_lrt, c(ip_mean_mle, 15))$root, error = function(e) NA))

    ip_sd_mle <- params_mle["ip_sd"]
    ip_sd_lrt <- function(ip_sd_now) ml - get.likelihood(c(ip_sd = ip_sd_now), data, params_mle) - qchisq(0.95, 1) / 2
    suppressWarnings(ip_sd_lower <- tryCatch(uniroot(ip_sd_lrt, c(1.5, ip_sd_mle))$root, error = function(e) NA))
    suppressWarnings(ip_sd_upper <- tryCatch(uniroot(ip_sd_lrt, c(ip_sd_mle, 5))$root, error = function(e) NA))

    data.frame(pi_mle = pi_mle,
               rho_mle = rho_mle,
               r_mle = r_mle,
               ip_mean_mle = ip_mean_mle,
               ip_sd_mle = ip_sd_mle,
               r_lower = r_lower,
               r_upper = r_upper,
               ip_mean_lower = ip_mean_lower,
               ip_mean_upper = ip_mean_upper,
               ip_sd_lower = ip_sd_lower,
               ip_sd_upper = ip_sd_upper)

}
