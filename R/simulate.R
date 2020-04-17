#' Simulate case information from the generative BETS model
#'
#' @keywords internal
#'
simulate.case <- function(n = 1e7,
                          params = c(pi = 0.1,
                                     lambda_w = 0.0005,
                                     lambda_v = 0.001,
                                     kappa = 4e-8,
                                     r = 0.2,
                                     alpha = 9,
                                     beta = 1.5,
                                     L = 53.5
                                     )) {

    param.names <- c("pi", "lambda_w", "lambda_v", "kappa", "r", "alpha", "beta", "L")
    stopifnot(all(param.names %in% names(params)))

    pi <- params["pi"]
    lambda_w <- params["lambda_w"]
    lambda_v <- params["lambda_v"]
    kappa <- params["kappa"]
    r <- params["r"]
    alpha <- params["alpha"]
    beta <- params["beta"]
    L <- params["L"]

    stopifnot(kappa < r * exp(- r * L))
    stopifnot(lambda_w * L < 1)
    stopifnot(lambda_v * L < 1)
    stopifnot(pi < 1)

    ## Begin of stay in Wuhan
    B <- sample(0:1, n, replace = TRUE, prob = c(1 - pi, pi))
    B[B > 0] <- runif(sum(B > 0), 0, L)

    ## End of stay in Wuhan
    E <- rep(NA, n)
    ind_wuhan <- (B == 0)
    E[ind_wuhan] <- sample(c(0, Inf), sum(ind_wuhan), replace = TRUE, prob = c(lambda_w * L, 1 - lambda_w * L))
    ind_visitor <- (B > 0)
    E[ind_visitor] <- sample(c(0, Inf), sum(ind_visitor), replace = TRUE, prob = c(lambda_v * L, 1 - lambda_v * L))
    ind_left <- (E == 0)
    E[ind_left] <- runif(sum(ind_left)) * L
    E[E < B] <- Inf

    ## Transmission
    FT <- runif(n)
    T <- log(FT * r / kappa) / r
    T[T < B | T > E] <- Inf

    ## Symptom
    S <- rep(Inf, n)
    ind_transmitted <- (T < Inf)
    S[ind_transmitted] <- T[ind_transmitted] + rgamma(sum(ind_transmitted), alpha, beta)

    ind <- (T < E) & (E < L) & (T < S) & (S < Inf)

    data.frame(B = B[ind],
               E = E[ind],
               T = T[ind],
               S = S[ind])

}
