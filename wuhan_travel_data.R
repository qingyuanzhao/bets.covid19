## Source: qianxi.baidu.com

date <- paste0(1:31, "-Jan")
arrive_rate <- c(2.85, 3.09, 4.22, 4.45, 5.08, 4.31, 4.25, 4.47, 4.81, 4.60, 4.64, 4.37, 4.83, 4.08, 4.06, 4.00, 4.40, 4.23, 4.15, 4.18, 4.24, 2.90, 1.75, 0.88, 0.63, 0.51, 0.42, 0.41, 0.37, 0.35, 0.33)
leave_rate <- c(3.46, 3.52, 5.52, 6.10, 5.32, 5.60, 6.41, 7.34, 8.14, 6.62, 7.56, 6.22, 5.76, 5.46, 5.91, 6.00, 6.44, 7.71, 7.41, 8.31, 10.74, 11.84, 11.14, 3.89, 1.30, 0.66, 0.43, 0.32, 0.26, 0.24, 0.24)
wuhan_travel <- data.frame(date = date, arrive_rate = arrive_rate, leave_rate = leave_rate)

## Source: https://www.jiqizhixin.com/articles/2020-01-27-2

taiwan_air <- 3696 + 2698 + 1121
japan_air <- 9080 + 6272 + 2656
singapore_air <- 10680
korea_air <- 6430
hongkong_air <- 7078
macau_air <- 6145

## Source: Wuhan Major (Zhou Xianwang)

WP_Jan26 <- 9000000
WP_Jan10 <- 14000000
