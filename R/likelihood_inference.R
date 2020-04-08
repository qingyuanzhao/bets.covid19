#' Approximate profile likelihood
#'
get.likelihood <- function(params, data, params_init = NULL, L = 54) {

    B <- data$B
    E <- data$E
    S <- data$S

    param.names <- c("rho", "r", "ip_mean", "ip_q95")

    if (all(param.names %in% names(params))) { # likelihood

        r <- params["r"]
        rho <- params["rho"]

        ## alpha <- (params["ip_mean"] / params["ip_sd"])^2
        ## beta <- params["ip_mean"] / params["ip_sd"]^2

        ## Maps (mean,q95) to (alpha,beta)
        param.map <- function(alpha, beta) {
            (alpha / beta - params["ip_mean"])^2 + (qgamma(0.95, alpha, beta) - params["ip_q95"])^2
        }
        ab <- optim(c(9, 1.5), function(ab) param.map(ab[1], ab[2]))$par
        alpha <- ab[1]
        beta <- ab[2]

        n <- length(B)

        ll <- 2 * n * log(r) + n * alpha * log(beta / (beta + r))
        ll <- ll - n * log(1 + rho * (1 - 2/r/L))
        for (i in 1:n) {
            if (B[i] > 0) {
                ll <- ll + log(rho / L)
            }
            ll <- ll + r * (S[i] - L)
            ll <- ll + log(pgamma(S[i] - B[i], alpha, beta + r) - pgamma(max(0, S[i] - E[i]), alpha, beta + r))
        }

        return(ll)
    } else { # profile likelihood

        params_init <- params_init[setdiff(param.names, names(params))]

        if (!all(param.names %in% names(c(params, params_init)))) {
            stop("To compute profile likelihood, params and params_init should together contain the following entries: r, ip_mean, ip_q95.")
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

    params_init <- c(rho = 1,
                     r = 0.2,
                     ip_mean = 5,
                     ip_q95 = 12)

    suppressWarnings(tmp <- optim(params_init,
                                  function(params) get.likelihood(params, data),
                                  control=list(fnscale=-1)))
    params_mle <- tmp$par
    ml <- tmp$value

    rho_mle <- params_mle["rho"]
    rho_lrt <- function(rho_now) ml - get.likelihood(c(rho = rho_now), data, params_mle) - qchisq(0.95, 1) / 2
    suppressWarnings(rho_lower <- tryCatch(uniroot(rho_lrt, c(rho_mle*0.2, rho_mle))$root, error = function(e) NA))
    suppressWarnings(rho_upper <- tryCatch(uniroot(rho_lrt, c(rho_mle, rho_mle*5))$root, error = function(e) NA))

    r_mle <- params_mle["r"]
    r_lrt <- function(r_now) ml - get.likelihood(c(r = r_now), data, params_mle) - qchisq(0.95, 1) / 2
    suppressWarnings(r_lower <- tryCatch(uniroot(r_lrt, c(0.1, r_mle))$root, error = function(e) NA))
    suppressWarnings(r_upper <- tryCatch(uniroot(r_lrt, c(r_mle, 0.6))$root, error = function(e) NA))

    ip_mean_mle <- params_mle["ip_mean"]
    ip_mean_lrt <- function(ip_mean_now) ml - get.likelihood(c(ip_mean = ip_mean_now), data, params_mle) - qchisq(0.95, 1) / 2
    suppressWarnings(ip_mean_lower <- tryCatch(uniroot(ip_mean_lrt, c(ip_mean_mle*0.5, ip_mean_mle))$root, error = function(e) NA))
    suppressWarnings(ip_mean_upper <- tryCatch(uniroot(ip_mean_lrt, c(ip_mean_mle, ip_mean_mle*1.5))$root, error = function(e) NA))

    ip_q95_mle <- params_mle["ip_q95"]
    ip_q95_lrt <- function(ip_q95_now) ml - get.likelihood(c(ip_q95 = ip_q95_now), data, params_mle) - qchisq(0.95, 1) / 2
    suppressWarnings(ip_q95_lower <- tryCatch(uniroot(ip_q95_lrt, c(ip_q95_mle*0.75, ip_q95_mle))$root, error = function(e) NA))
    suppressWarnings(ip_q95_upper <- tryCatch(uniroot(ip_q95_lrt, c(ip_q95_mle, ip_q95_mle*1.5))$root, error = function(e) NA))

    data.frame(rho_mle = rho_mle,
               rho_lower = rho_lower,
               rho_upper = rho_upper,
               r_mle = r_mle,
               r_lower = r_lower,
               r_upper = r_upper,
               ip_mean_mle = ip_mean_mle,
               ip_mean_lower = ip_mean_lower,
               ip_mean_upper = ip_mean_upper,
               ip_q95_mle = ip_q95_mle,
               ip_q95_lower = ip_q95_lower,
               ip_q95_upper = ip_q95_upper)

}
