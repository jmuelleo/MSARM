# Refits MSARM to audit variance-floor hits against the original fit scores.
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Supply batch 1-12 and an output CSV.")
source("benchmarks/stress-regimes.R")
fits <- read.csv("benchmarks/stress300-results.csv")
fits <- fits[fits$package == "MSARM", ]
batch <- as.integer(args[1L])
indices <- seq.int(batch, nrow(fits), by = 12L)
rows <- lapply(indices, function(i) {
  row <- fits[i, ]
  s <- scenarios[[row$scenario]]
  data <- simulate(s, 1000L + row$replication)
  set.seed(2000L + row$replication)
  f <- env$MSARM.fit(data$y, data$order, s$k, s$sw, R_value = 5)
  stopifnot(abs(f$LV - row$loglik) < 1e-7)
  variance <- f$alpha[, data$order + 2L]
  minimum <- which.min(variance)
  sample_var <- var(data$y)
  data.frame(scenario = row$scenario, replication = row$replication,
             min_variance = variance[minimum],
             floor_ratio = variance[minimum] / (sample_var * 1e-10),
             effective_regime_count = sum(f$zeta_t_T[minimum, ]),
             wrong = row$wrong, failure = row$failure)
})
write.csv(do.call(rbind, rows), args[2L], row.names = FALSE)
