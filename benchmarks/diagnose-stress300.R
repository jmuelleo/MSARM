# Refits named scenario/replication pairs and inspects all five MSwM starts.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply an output CSV and scenario/replication pairs.")
source("benchmarks/stress-regimes.R")
describe <- function(P, beta, s, data) {
  k <- s$k
  K <- data$order
  probabilities <- env$zeta_YT_function(beta, P, data$y, k, K)
  labels <- classify(probabilities[, -seq_len(K), drop = FALSE],
                     data$state[-seq_len(K)], k)
  c(list(loglik = score(data$y, K, P, beta),
         intercepts = paste(signif(beta[, 1L], 5), collapse = ":"),
         phi1 = paste(signif(beta[, 2L], 5), collapse = ":"),
         variances = paste(signif(beta[, K + 2L], 5), collapse = ":"),
         pstay = paste(signif(diag(P), 5), collapse = ":")),
    labels[c("wrong", "balanced_wrong", "true_counts", "predicted_counts")])
}
refine <- function(P, beta, s, data) {
  design <- env$msarm_design(data$y, data$order)
  for (i in seq_len(250L)) {
    filter <- env$msarm_filter(beta, P, data$y, s$k, data$order, design)
    smooth <- env$msarm_smooth(filter, P, data$order)
    beta <- env$msarm_update(design$y, design$X,
      smooth[, -seq_len(data$order), drop = FALSE], beta, s$sw,
      var(data$y) * 1e-10)
    P <- env$msarm_transition(P, filter, smooth, data$order, 1L)
  }
  describe(P, beta, s, data)
}
rows <- list()
for (item in args[-1L]) {
  id <- strsplit(item, "/", fixed = TRUE)[[1L]]
  s <- scenarios[[id[1L]]]
  i <- as.integer(id[2L])
  data <- simulate(s, 1000L + i)
  output <- function(method, start, fit = NULL, error = "") {
    if (is.null(fit)) {
      fit <- as.list(rep(NA, 9L))
      names(fit) <- c("loglik", "intercepts", "phi1", "variances", "pstay",
                      "wrong", "balanced_wrong", "true_counts", "predicted_counts")
    }
    as.data.frame(c(list(scenario = id[1L], replication = i,
                        method = method, start = start, error = error), fit))
  }
  if (is.null(s$innovation)) {
    rows[[length(rows) + 1L]] <- output("generating", NA,
      describe(data$P, s$beta, s, data))
  }
  set.seed(2000L + i)
  f <- env$MSARM.fit(data$y, data$order, s$k, s$sw, R_value = 5)
  rows[[length(rows) + 1L]] <- output("MSARM", NA,
    describe(f$P, f$alpha, s, data))
  set.seed(2000L + i)
  base <- lm(y ~ 1, data = data.frame(y = data$y))
  best <- NULL
  best_ll <- -Inf
  for (start in seq_len(5L)) {
    result <- tryCatch({
      invisible(capture.output(model <- suppressWarnings(MSwM::msmFit(
        base, k = s$k, sw = s$sw, p = data$order,
        control = list(parallelization = FALSE, maxiter = 250)))))
      beta <- cbind(as.matrix(model@Coef), model@std^2)
      metrics <- describe(model@transMat, beta, s, data)
      if (is.finite(metrics$loglik) && metrics$loglik > best_ll) {
        best <- list(P = model@transMat, beta = beta)
        best_ll <- metrics$loglik
      }
      output("MSwM", start, metrics)
    }, error = function(e) output("MSwM", start, error = conditionMessage(e)))
    rows[[length(rows) + 1L]] <- result
  }
  if (!is.null(best)) {
    rows[[length(rows) + 1L]] <- tryCatch(
      output("MSARM from best MSwM", NA, refine(best$P, best$beta, s, data)),
      error = function(e) output("MSARM from best MSwM", NA,
                                 error = conditionMessage(e)))
  }
  write.csv(do.call(rbind, rows), args[1L], row.names = FALSE)
  message("Diagnosed ", item)
}
