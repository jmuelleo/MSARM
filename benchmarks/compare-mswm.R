# Run from the repository root with Rscript benchmarks/compare-mswm.R 16.
suppressPackageStartupMessages(library(MSwM))
env <- new.env()
invisible(lapply(list.files("R", full.names = TRUE), sys.source, envir = env))
args <- commandArgs(trailingOnly = TRUE)
replications <- if (length(args)) as.integer(args[1L]) else 16L

scenarios <- list(
  intercept = list(k = 2L, n = 250L, sw = c(TRUE, FALSE, FALSE),
                   beta = rbind(c(-1.5, .3, 1), c(1.5, .3, 1))),
  full = list(k = 2L, n = 400L, sw = c(TRUE, TRUE, TRUE),
              beta = rbind(c(-1.5, .65, .5), c(1.5, -.45, 1.4))),
  variance = list(k = 2L, n = 350L, sw = c(FALSE, FALSE, TRUE),
                  beta = rbind(c(0, .4, .2), c(0, .4, 2))),
  mixed_ar2 = list(k = 2L, n = 350L, sw = c(TRUE, FALSE, TRUE, TRUE),
                   beta = rbind(c(-1.2, .3, -.2, .4), c(1.2, .3, .25, 1))),
  weak_short = list(k = 2L, n = 80L, sw = c(TRUE, TRUE, FALSE),
                     beta = rbind(c(-.35, .3, 1), c(.35, .5, 1))),
  three_ar2 = list(k = 3L, n = 500L, sw = c(TRUE, TRUE, FALSE, TRUE),
                    beta = rbind(c(-2, .4, -.1, .4), c(0, -.5, -.1, 1),
                                 c(2, .55, -.1, 1.4)))
)
simulate <- function(s, seed) {
  set.seed(seed)
  P <- matrix(.08 / (s$k - 1L), s$k, s$k)
  diag(P) <- .92
  order <- ncol(s$beta) - 2L
  y <- numeric(s$n + 100L)
  state <- rep(1L, length(y))
  for (i in seq.int(order + 1L, length(y))) {
    j <- state[i] <- sample.int(s$k, 1L, prob = P[, state[i - 1L]])
    y[i] <- sum(c(1, y[i - seq_len(order)]) * s$beta[j, -ncol(s$beta)]) +
      rnorm(1, sd = sqrt(s$beta[j, ncol(s$beta)]))
  }
  list(y = tail(y, s$n), state = tail(state, s$n), P = P, order = order)
}
# An independent density-scale recursion scores both fits with the same
# initial distribution, lag alignment, and treatment of the first observation.
score <- function(y, order, P, beta) {
  lagged <- embed(y, order + 1L)
  X <- cbind(1, lagged[, -1L, drop = FALSE])
  means <- X %*% t(beta[, -ncol(beta), drop = FALSE])
  density <- matrix(dnorm(rep(lagged[, 1L], nrow(beta)), means,
                          rep(sqrt(beta[, ncol(beta)]), each = nrow(X))),
                     nrow(X))
  pred <- rep(1 / nrow(P), nrow(P))
  value <- 0
  for (i in seq_len(nrow(X))) {
    w <- pred * density[i, ]
    z <- sum(w)
    if (!is.finite(z) || z <= 0) return(-Inf)
    value <- value + log(z)
    pred <- drop(P %*% (w / z))
  }
  value
}
fit_one <- function(package, data, s, seed) {
  set.seed(seed)
  y <- data$y
  start <- proc.time()[["elapsed"]]
  result <- tryCatch({
    candidates <- if (package == "MSARM") {
      fit <- env$MSARM.fit(y, data$order, s$k, s$sw, R_value = 5)
      list(list(P = fit$P, beta = fit$alpha))
    } else {
      base <- lm(y ~ 1, data = data.frame(y))
      lapply(seq_len(5L), function(i) tryCatch({
        invisible(capture.output(fit <- suppressWarnings(MSwM::msmFit(
          base, k = s$k, sw = s$sw, p = data$order,
          control = list(parallelization = FALSE, maxiter = 250)))))
        list(P = fit@transMat, beta = cbind(as.matrix(fit@Coef), fit@std^2))
      }, error = identity))
    }
    values <- vapply(candidates, function(fit) {
      if (inherits(fit, "error")) return(-Inf)
      score(y, data$order, fit$P, fit$beta)
    }, 0)
    if (!any(is.finite(values))) stop("No finite score from five starts")
    best <- candidates[[which.max(values)]]
    smooth <- env$zeta_YT_function(best$beta, best$P, y, s$k, data$order)
    guess <- max.col(t(smooth[, -seq_len(data$order)]), ties.method = "first")
    truth <- data$state[-seq_len(data$order)]
    permutations <- if (s$k == 2L) rbind(1:2, 2:1) else
      rbind(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
    error_rate <- min(apply(permutations, 1L, function(p) mean(p[guess] != truth)))
    list(loglik = max(values), misclassified = error_rate, error = "")
  }, error = function(e) list(loglik = NA_real_, misclassified = NA_real_,
                              error = conditionMessage(e)))
  c(list(package = package, seconds = proc.time()[["elapsed"]] - start), result)
}
selected <- if (length(args) >= 2L) args[2L] else names(scenarios)
if (!all(selected %in% names(scenarios))) stop("Unknown scenario")
output <- if (length(args) >= 3L) args[3L] else "benchmarks/results.csv"
offset <- if (length(args) >= 4L) as.integer(args[4L]) else 0L
rows <- list()
for (name in selected) {
  s <- scenarios[[name]]
  for (i in offset + seq_len(replications)) {
    data <- simulate(s, 1000L + i)
    truth <- score(data$y, data$order, data$P, s$beta)
    for (package in c("MSARM", "MSwM")) {
      fit <- fit_one(package, data, s, 2000L + i)
      rows[[length(rows) + 1L]] <- as.data.frame(
        c(list(scenario = name, replication = i, truth = truth), fit))
    }
    message("Completed ", name, " replication ", i)
  }
  message("Completed ", name)
}
results <- do.call(rbind, rows)
write.csv(results, output, row.names = FALSE)
print(aggregate(cbind(seconds, misclassified) ~ package, results, median))
print(with(results, table(package, failed = !is.finite(loglik))))
