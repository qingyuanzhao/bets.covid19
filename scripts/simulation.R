devtools::load_all()

params_true <- c(pi = 0.1,
                 lambda_w = 0.0005,
                 lambda_v = 0.001,
                 kappa = 4e-8,
                 r = 0.2,
                 alpha = 9,
                 beta = 1.5,
                 L = 54)


res <- parallel::mclapply(1:30, function(sim) {
    print(sim)
    data <- simulate.case(params = params_true)
    likelihood.inference(data)
    }, mc.cores = 3)

res <- do.call(rbind, res)

1 - mean(0.25 < res$r_lower | 0.25 > res$r_upper, na.rm = TRUE)
1 - mean(12 / 1.5 < res$ip_mean_lower | 12/1.5 > res$ip_mean_upper, na.rm = TRUE)
1 - mean(sqrt(12)/1.5 < res$ip_sd_lower | sqrt(12)/1.5 > res$ip_sd_upper, na.rm = TRUE)
