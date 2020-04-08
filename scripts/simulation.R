devtools::load_all()

params_true <- c(pi = 0.2,
                 lambda_w = 0.002,
                 lambda_v = 0.004,
                 kappa = 1e-10,
                 r = 0.3,
                 alpha = 2,
                 beta = 0.47,
                 L = 54) # nu = 1

## check inclusion probability is approximately correctly
tmp <- as.list(params_true)
attach(tmp)
(1-pi) * lambda_w * kappa / r^2 * exp(r * L)
pi * lambda_v * kappa / r^2 * exp(r * L) * (1 - 2 / r / L)
detach(tmp)

inclusion_w <- replicate(10 , {data <- simulate.case(params = params_true); sum(data$B == 0) / 1e7})
mean(inclusion_w)
inclusion_v <- replicate(10 , {data <- simulate.case(params = params_true); sum(data$B > 0) / 1e7})
mean(inclusion_v)

## Simulations

res <- parallel::mclapply(1:30, function(sim) {
    print(sim)
    data <- simulate.case(params = params_true)
    print(res_one <- likelihood.inference(data))
    res_one
    }, mc.cores = 3)

res <- do.call(rbind, res)

1 - mean(params_true["r"] < res$r_lower | params_true["r"] > res$r_upper, na.rm = TRUE)
## 1 - mean(params_true$alpha / params_true$beta < res$ip_mean_lower | params_true$alpha / params_true$beta > res$ip_mean_upper, na.rm = TRUE)
