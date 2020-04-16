#' Approximate profile likelihood
#'
#' @param params A vector of parameters (names: rho, r, ip_q50, ip_q95)
#' @param data A data frame with three columns: B, E, S
#' @param params_init Initial parameters for computing the profile likelihood
#' @param L Day of travel quarantine
#'
#' @return When \code{params} contains all the parameters (rho, r, ip_q50, ip_q95), returns the approximate log-likelihood of \code{data}. When \code{params} contains some but not all the parameters, returns the profile log-likelihood.
#'
get.likelihood.unconditional <- function(params, data, params_init = NULL, L = 54) {

    B <- data$B
    E <- data$E
    S <- data$S

    param.names <- c("rho", "r", "ip_q50", "ip_q95")

    if (all(param.names %in% names(params))) { # likelihood

        r <- params["r"]
        rho <- params["rho"]

        ## alpha <- (params["ip_q50"] / params["ip_sd"])^2
        ## beta <- params["ip_q50"] / params["ip_sd"]^2

        ## Maps (mean,q95) to (alpha,beta)
        param.map <- function(alpha, beta) {
            (qgamma(0.5, alpha, beta) - params["ip_q50"])^2 + (qgamma(0.95, alpha, beta) - params["ip_q95"])^2
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
            stop("To compute profile likelihood, params and params_init should together contain the following entries: r, ip_q50, ip_q95.")
        }

        optim(params_init,
              function(params_now) get.likelihood.unconditional(c(params, params_now), data),
              control=list(fnscale=-1))$value

    }

}


#' (Profile) Conditional likelihood given B and E
#'
#' @param params A vector of parameters (names: r, ip_q50, ip_q95)
#' @param data A data frame with three columns: B, E, S
#' @param M Right truncation for symptom onset
#' @param r Parameter for epidemic growth (overrides \code{{params}})
#' @param params_init Initial parameters for computing the profile likelihood
#'
#' @return Conditional log-likelihood.
#'
#'
get.likelihood.conditional <- function(params, data, M = Inf, r = NULL, params_init = NULL) {

    if (!is.null(r)) {
        params["r"] <- r
    }

    param.names <- c("r", "ip_q50", "ip_q95")

    if (all(param.names %in% names(params))) { # likelihood

        r <- params["r"]

        param.map <- function(alpha, beta) {
            (qgamma(0.5, alpha, beta) - params["ip_q50"])^2 + (qgamma(0.95, alpha, beta) - params["ip_q95"])^2
        }
        ab <- optim(c(9, 1.5), function(ab) param.map(ab[1], ab[2]))$par
        alpha <- ab[1]
        beta <- ab[2]

        if (M == Inf) {

            B <- data$B
            E <- data$E
            S <- data$S
            n <- length(B)

            if (r == 0) {
                ll <- 0
                for (i in 1:n) {
                    ll <- ll - log(E[i] - B[i])
                    ll <- ll + log(pgamma(S[i] - B[i], alpha, beta) - pgamma(max(0, S[i] - E[i]), alpha, beta))
                }
            } else {
                ll <- n * log(r) + n * alpha * log(beta / (beta + r))
                for (i in 1:n) {
                    ll <- ll + r * S[i] + log(pgamma(S[i] - B[i], alpha, beta + r) - pgamma(max(0, S[i] - E[i]), alpha, beta + r))
                    ll <- ll - log(exp(r * E[i]) - exp(r * B[i]))
                }
            }
        } else { # Adjust for ascertainment
            data <- subset(data, S <= M)
            B <- data$B
            E <- data$E
            S <- data$S
            n <- length(B)

            if (r == 0) {
                myfun <- function(x) {
                    x * pgamma(x, alpha, beta) - (akoha / beta) * pgamma(x, alpha + 1, beta)
                }
                ll <- 0
                for (i in 1:n) {
                    ll <- ll + log(pgamma(S[i] - B[i], alpha, beta) - pgamma(max(0, S[i] - E[i]), alpha, beta))
                    ll <- ll - log(myfun(M - B[i]) - myfun(max(0, M - E[i])))
                }
            } else {
                myfun <- function(x) {
                    (beta / (beta + r))^alpha * pgamma(x, alpha, beta + r) - exp(- r * x) * pgamma(x, alpha, beta)
                }
                ll <- n * log(r) + n * alpha * log(beta / (beta + r))
                for (i in 1:n) {
                    ll <- ll + r * (S[i] - M) + log(pgamma(S[i] - B[i], alpha, beta + r) - pgamma(max(0, S[i] - E[i]), alpha, beta + r))
                    ll <- ll - log(myfun(M - B[i]) - myfun(max(0, M - E[i])))
                }
            }
        }
    } else { # profile likelihood

        params_init <- params_init[setdiff(param.names, names(params))]

        if (!all(param.names %in% names(c(params, params_init)))) {
            stop("To compute profile likelihood, params and params_init should together contain the following entries: r, ip_q50, ip_q95.")
        }

        ll <- optim(params_init,
                    function(params_now) get.likelihood.conditional(c(params, params_now), data, M = M, r = r),
                    control=list(fnscale=-1))$value

    }

    return(ll)

}

#' Likelihood function
#'
#'
get.likelihood <- function(params, data, likelihood = c("conditional", "unconditional"), M = Inf, r = NULL, params_init = NULL, L = 54) {

    likelihood <- match.arg(likelihood, c("conditional", "unconditional"))
    switch(likelihood,
           conditional = get.likelihood.conditional(params, data, M = M, r = r, params_init = params_init),
           unconditional = get.likelihood.unconditional(params, data, params_init = params_init, L = L))

}

#' Likelihood inference
#'
#' @param data A data.frame with three columns: B, E, S.
#' @param likelihood Conditional on B and E?
#' @param ci How to compute the confidence interval?
#' @param M Right truncation for symptom onset (only available for conditional likelihood)
#' @param r Parameter for epidemic growth (overrides \code{{params}, only available for conditional likelihood})
#' @param level Level of the confidence interval (default 0.95).
#' @param bootstrap Number of bootstrap resamples.
#' @param mc.cores Number of cores used for computing the bootstrap confidence interval.
#'
#' @details The confidence interval is either not computed (\code{"point"}), or computed by inverting the likelihood ratio test (\code{"lrt"}) or basic bootstrap (\code{"bootstrap"})
#'
#' @import rootSolve
#' @importFrom parallel mclapply
#'
#' @export
#'
likelihood.inference <- function(data, likelihood = c("conditional", "unconditional"), ci = c("lrt", "point", "bootstrap"),M = Inf, r = NULL, L = 54, level = 0.95, bootstrap = 1000, mc.cores = 1) {

    likelihood <- match.arg(likelihood, c("conditional", "unconditional"))
    ci <- match.arg(ci, c("lrt", "point", "bootstrap"))

    if (likelihood == "unconditional") {
        params_init <- c(rho = 1,
                         r = 0.2,
                         ip_q50 = 5,
                         ip_q95 = 12)
    } else if (is.null(r)) {
        params_init <- c(r = 0.2,
                         ip_q50 = 5,
                         ip_q95 = 12)
    } else {
        params_init <- c(ip_q50 = 5,
                         ip_q95 = 12)
    }

    suppressWarnings(tmp <- optim(params_init,
                                  function(params) get.likelihood(params, data, likelihood = likelihood, M = M, r = r, L = L),
                                  control=list(fnscale=-1)))
    params_mle <- tmp$par
    ml <- tmp$value

    if (likelihood == "unconditional") {
        rho_mle <- params_mle["rho"]
    } else {
        rho_mle <- NA
    }
    if (is.null(r)) {
        r_mle <- params_mle["r"]
    } else {
        r_mle <- r
    }
    ip_q50_mle <- params_mle["ip_q50"]
    ip_q95_mle <- params_mle["ip_q95"]

    if (ci == "point") {
        res <- data.frame(sample_size = nrow(data),
                          rho_mle = rho_mle,
                          r_mle = r_mle,
                          ip_q50_mle = ip_q50_mle,
                          ip_q95_mle = ip_q95_mle)
        rownames(res) <- NULL

        return(res)

    } else if (ci == "lrt") {
        if (likelihood == "unconditional") {
            rho_lrt <- function(rho_now) ml - get.likelihood(c(rho = rho_now), data, likelihood = likelihood, M = M, r = r, L = L, params_init = params_mle) - qchisq(level, 1) / 2
            suppressWarnings(rho_lower <- tryCatch(uniroot(rho_lrt, c(rho_mle*0.2, rho_mle))$root, error = function(e) NA))
            suppressWarnings(rho_upper <- tryCatch(uniroot(rho_lrt, c(rho_mle, rho_mle*5))$root, error = function(e) NA))
        } else {
            rho_lower <- rho_upper <- NA
        }

        if (is.null(r)) {
            r_lrt <- function(r_now) ml - get.likelihood(c(r = r_now), data, likelihood = likelihood, M = M, r = r, L = L, params_init = params_mle) - qchisq(level, 1) / 2
            suppressWarnings(r_lower <- tryCatch(uniroot(r_lrt, c(0.1, r_mle))$root, error = function(e) NA))
            suppressWarnings(r_upper <- tryCatch(uniroot(r_lrt, c(r_mle, 0.8))$root, error = function(e) NA))
        } else {
            r_lower <- r_upper <- NA
        }

        ip_q95_lrt <- function(ip_q95_now) ml - get.likelihood(c(ip_q95 = ip_q95_now), data, likelihood = likelihood, M = M, r = r, L = L, params_init = params_mle) - qchisq(level, 1) / 2
        suppressWarnings(ip_q95_lower <- tryCatch(uniroot(ip_q95_lrt, c(ip_q95_mle*0.75, ip_q95_mle))$root, error = function(e) NA))
        suppressWarnings(ip_q95_upper <- tryCatch(uniroot(ip_q95_lrt, c(ip_q95_mle, ip_q95_mle*1.5))$root, error = function(e) NA))

        ip_q50_lrt <- function(ip_q50_now) ml - get.likelihood(c(ip_q50 = ip_q50_now), data, likelihood = likelihood, M = M, r = r, L = L, params_init = params_mle) - qchisq(level, 1) / 2
        suppressWarnings(ip_q50_lower <- tryCatch(uniroot(ip_q50_lrt, c(ip_q50_mle*0.5, ip_q50_mle))$root, error = function(e) NA))
        suppressWarnings(ip_q50_upper <- tryCatch(uniroot(ip_q50_lrt, c(ip_q50_mle, ip_q50_mle*1.5))$root, error = function(e) NA))

    } else {

        res_bootstrap <- parallel::mclapply(1:bootstrap, function(iter) {
            if (iter %% 50 == 0) {
                print(paste("Bootstrap resample:", iter))
            }
            likelihood.inference(data[sample(1:nrow(data), replace = TRUE), ], likelihood, ci = "point", M = M, r = r, L = L)
        }, mc.cores = mc.cores)

        res_bootstrap <- do.call(rbind, res_bootstrap)

        rho_lower <- 2 * rho_mle - quantile(res_bootstrap$rho_mle, 1 - (1 - level) / 2, na.rm = TRUE)
        rho_upper <- 2 * rho_mle - quantile(res_bootstrap$rho_mle, (1 - level) / 2, na.rm = TRUE)

        r_lower <- 2 * r_mle - quantile(res_bootstrap$r_mle, 1 - (1 - level) / 2, na.rm = TRUE)
        r_upper <- 2 * r_mle - quantile(res_bootstrap$r_mle, (1 - level) / 2, na.rm = TRUE)

        ip_q50_lower <- 2 * ip_q50_mle - quantile(res_bootstrap$ip_q50_mle, 1 - (1 - level) / 2, na.rm = TRUE)
        ip_q50_upper <- 2 * ip_q50_mle - quantile(res_bootstrap$ip_q50_mle, (1 - level) / 2, na.rm = TRUE)

        ip_q95_lower <- 2 * ip_q95_mle - quantile(res_bootstrap$ip_q95_mle, 1 - (1 - level) / 2, na.rm = TRUE)
        ip_q95_upper <- 2 * ip_q95_mle - quantile(res_bootstrap$ip_q95_mle, (1 - level) / 2, na.rm = TRUE)

    }

    res <- data.frame(sample_size = nrow(data),
                      rho_mle = rho_mle,
                      rho_lower = rho_lower,
                      rho_upper = rho_upper,
                      r_mle = r_mle,
                      r_lower = r_lower,
                      r_upper = r_upper,
                      ip_q50_mle = ip_q50_mle,
                      ip_q50_lower = ip_q50_lower,
                      ip_q50_upper = ip_q50_upper,
                      ip_q95_mle = ip_q95_mle,
                      ip_q95_lower = ip_q95_lower,
                      ip_q95_upper = ip_q95_upper)
    rownames(res) <- NULL

    return(res)

}
