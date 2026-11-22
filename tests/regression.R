library(MSARM)

P <- matrix(c(0.85, 0.15, 0.25, 0.75), 2L)
alpha <- rbind(c(-1, 0.6, 0.5), c(1, -0.4, 1.5))
y <- c(-0.4, 0.2, 1.5, -0.2)
paths <- as.matrix(expand.grid(rep(list(1:2), length(y) - 1L)))
weights <- apply(paths, 1L, function(s) {
  densities <- dnorm(y[-1L], alpha[s, 1L] + alpha[s, 2L] * head(y, -1L),
                     sqrt(alpha[s, 3L]))
  0.5 * prod(densities) * prod(P[cbind(tail(s, -1L), head(s, -1L))])
})
smooth <- sapply(seq_len(ncol(paths)), function(i) {
  vapply(1:2, function(j) sum(weights[paths[, i] == j]) / sum(weights), 0)
})
stopifnot(isTRUE(all.equal(Log_Likelihood_function(alpha, P, y, 2, 1),
                          log(sum(weights)))),
          isTRUE(all.equal(zeta_YT_function(alpha, P, y, 2, 1)[, -1L],
                          smooth)))
transitions <- Reduce(`+`, lapply(seq_len(nrow(paths)), function(i) {
  s <- paths[i, ]
  matrix(tabulate(tail(s, -1L) + 2L * (head(s, -1L) - 1L), 4L), 2L) *
    weights[i]
}))
stopifnot(isTRUE(all.equal(P_ij_next_function(alpha, P, y, 2, 1, 1),
                          sweep(transitions, 2L, colSums(transitions), "/"))))

fixture <- list(P = P, alpha = alpha, zeta_t_T = smooth, Y_T = y,
                residuals = c(-1, 0, 1))
expected <- numeric(3L)
for (initial in 1:2) {
  for (i in seq_len(nrow(paths))) {
    regimes <- paths[i, ]
    previous <- initial
    value <- tail(y, 1L)
    weight <- tail(t(smooth), 1L)[initial]
    forecasts <- numeric(3L)
    for (h in seq_len(3L)) {
      current <- regimes[h]
      weight <- weight * P[current, previous]
      value <- alpha[current, 1L] + alpha[current, 2L] * value
      forecasts[h] <- value
      previous <- current
    }
    expected <- expected + weight * forecasts
  }
}
set.seed(13)
forecast <- MSARM.predict(fixture, 3, boot = TRUE, L = 10000)
stopifnot(isTRUE(all.equal(forecast$fcast, expected)),
          max(abs(forecast$fcast_boot[, 1L] - expected)) < 0.08,
          all(abs(colSums(forecast$zeta_hat) - 1) < 1e-10))

set.seed(15)
y <- arima.sim(list(ar = c(0.5, -0.2)), n = 90)
patterns <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 4L)))
for (i in seq_len(nrow(patterns))) {
  sw <- patterns[i, ]
  set.seed(19)
  first <- MSARM.fit(y, 2, 3, sw, max_value = 3, R_value = 1)
  set.seed(19)
  next_fit <- MSARM.fit(y, 2, 3, sw, max_value = 4, R_value = 1)
  stopifnot(next_fit$LV >= first$LV - 1e-8,
            isTRUE(all.equal(next_fit$LV, Log_Likelihood_function(
              next_fit$alpha, next_fit$P, y, 3, 2))),
            max(abs(colSums(next_fit$P) - 1)) < 1e-10,
            max(abs(colSums(next_fit$zeta_t_T) - 1)) < 1e-10,
            all(apply(next_fit$alpha[, !sw, drop = FALSE], 2L,
                      function(x) all(x == x[1L]))))
}

set.seed(21)
fit <- MSARM.fit(y, 2, 2, c(TRUE, FALSE, TRUE, TRUE),
                 max_value = 15, R_value = 1)
set.seed(21)
scaled <- MSARM.fit(y * 1e7, 2, 2, c(TRUE, FALSE, TRUE, TRUE),
                    max_value = 15, R_value = 1)
stopifnot(isTRUE(all.equal(fit$P, scaled$P, tolerance = 1e-6)),
          isTRUE(all.equal(fit$LV, scaled$LV + (length(y) - 2) * log(1e7),
                           tolerance = 1e-6)))

# Zero predicted mass and a one-observation post-lag sample.
absorbing <- rbind(c(1, 1), c(0, 0))
single <- zeta_YT_function(alpha, absorbing, c(0.1, 0.2), 2, 1)
stopifnot(isTRUE(all.equal(sum(single[, 2L]), 1)))
zero <- zeta_YT_function(alpha, absorbing, c(0.1, 0.2, 1, 2), 2, 1)
stopifnot(all(is.finite(zero[, -1L])), all(zero[2L, -c(1L, 2L)] == 0))

# Short time series forecast bands and one-path bootstrap.
fixture$Y_T <- ts(fixture$Y_T, start = c(2000, 2), frequency = 4)
one <- MSARM.predict(fixture, boot = TRUE, levels = 0.9, L = 1)
stopifnot(identical(dim(one$fcast_boot), c(1L, 3L)),
          all(one$fcast_boot == one$fcast_boot[1L]))
grDevices::pdf(NULL)
MSARM.plot(one, conf = TRUE)
invisible(grDevices::dev.off())

# The three-regime RCM is bounded.
flat <- MSARM.fit(y, 2, 3, rep(FALSE, 4), max_value = 3, R_value = 1)
stopifnot(flat$RCM >= 0, flat$RCM <= 100 + 1e-10)
