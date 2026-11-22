# Run from the repository root with Rscript benchmarks/stress-regimes.R 16.
source("benchmarks/compare-mswm.R", local = TRUE)
make_case <- function(parent, n, P, beta = NULL) {
  s <- scenarios[[parent]]
  s$n <- n
  s$P <- P
  if (!is.null(beta)) s$beta <- beta
  s
}
scenarios$rare <- make_case("intercept", 250L,
                            rbind(c(.985, .30), c(.015, .70)))
scenarios$brief <- make_case("intercept", 250L,
                             rbind(c(.96, .72), c(.04, .28)))
scenarios$fast <- make_case("full", 250L,
                            rbind(c(.60, .40), c(.40, .60)))
scenarios$weak_long <- make_case("weak_short", 250L,
                                 rbind(c(.97, .03), c(.03, .97)))
scenarios$near_unit <- make_case("full", 250L,
                                 rbind(c(.94, .06), c(.06, .94)),
                                 rbind(c(-.2, .97, .6), c(.4, .8, 1.3)))
scenarios$rare_three <- make_case("three_ar2", 250L,
  rbind(c(.95, .05, .20), c(.045, .94, .25), c(.005, .01, .55)))

# Additional scenarios used in the 300-series study.
scenarios$clear <- make_case("intercept", 160L, rbind(c(.96,.04), c(.04,.96)),
                             rbind(c(-3,.2,.5), c(3,.2,.5)))
scenarios$rare_clear <- make_case("intercept", 180L, rbind(c(.98,.25), c(.02,.75)),
                                  rbind(c(-2,.2,.5), c(3,.2,.5)))
scenarios$rare_weak <- make_case("intercept", 160L, rbind(c(.98,.30), c(.02,.70)),
                                 rbind(c(-.2,.3,1), c(.6,.3,1)))
scenarios$overlap <- make_case("intercept", 150L, rbind(c(.90,.10), c(.10,.90)),
                               rbind(c(-.2,.3,1), c(.2,.3,1)))
scenarios$ar_only <- make_case("full", 180L, rbind(c(.94,.06), c(.06,.94)),
                               rbind(c(0,.8,1), c(0,-.7,1)))
scenarios$ar_only$sw <- c(FALSE, TRUE, FALSE)
scenarios$short_clear <- make_case("clear", 55L, rbind(c(.90,.10), c(.10,.90)))
scenarios$short_var <- make_case("variance", 70L, rbind(c(.92,.08), c(.08,.92)))
scenarios$quick_var <- make_case("variance", 160L, rbind(c(.56,.44), c(.44,.56)))
scenarios$clear_three <- make_case("three_ar2", 155L,
  rbind(c(.92,.04,.04), c(.04,.92,.04), c(.04,.04,.92)),
  rbind(c(-4,.25,.5), c(0,.25,.5), c(4,.25,.5)))
scenarios$clear_three$sw <- c(TRUE, FALSE, FALSE)
scenarios$near_absorb <- make_case("clear", 175L, rbind(c(.995,.005), c(.005,.995)))
scenarios$heavy_tail <- make_case("intercept", 180L, rbind(c(.94,.06), c(.06,.94)))
scenarios$heavy_tail$innovation <- "t3"
scenarios$outliers <- make_case("intercept", 180L, rbind(c(.94,.06), c(.06,.94)))
scenarios$outliers$innovation <- "outliers"
scenarios$long_weak <- make_case("weak_short", 500L, rbind(c(.96,.04), c(.04,.96)))

classify <- function(smooth, truth, k) {
  guess <- max.col(t(smooth), ties.method = "first")
  permutations <- if (k == 2L) rbind(1:2, 2:1) else
    rbind(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  actual <- tabulate(truth, k)
  matrices <- lapply(seq_len(nrow(permutations)), function(i) {
    aligned <- permutations[i, guess]
    matrix(tabulate(truth + k * (aligned - 1L), k * k), k)
  })
  correct <- vapply(matrices, function(m) sum(diag(m)), 0L)
  balanced <- vapply(matrices, function(m) {
    mean(diag(m)[actual > 0L] / actual[actual > 0L])
  }, 0)
  m <- matrices[[order(-correct, -balanced)[1L]]]
  predicted <- colSums(m)
  recall <- diag(m) / ifelse(actual == 0L, NA, actual)
  list(wrong = 1 - sum(diag(m)) / length(truth),
       balanced_wrong = 1 - mean(recall, na.rm = TRUE),
       min_recall = min(recall, na.rm = TRUE),
       true_states = sum(actual > 0L), predicted_states = sum(predicted > 0L),
       collapsed = sum(actual > 0L) > 1L && sum(predicted > 0L) == 1L,
       never_predicted = sum(actual > 0L & predicted == 0L),
       zero_recall = sum(recall == 0, na.rm = TRUE),
       true_counts = paste(actual, collapse = ":"),
       predicted_counts = paste(predicted, collapse = ":"),
       confusion = paste(c(t(m)), collapse = ":"))
}
# Exercise the low-error collapse case and arbitrary label reversal.
case <- classify(rbind(rep(1, 100), rep(0, 100)), c(rep(1L, 98), 2L, 2L), 2L)
stopifnot(case$collapsed, case$never_predicted == 1, case$zero_recall == 1,
          abs(case$wrong - .02) < 1e-12, case$balanced_wrong == .5)
case <- classify(rbind(c(0, 1), c(1, 0)), 1:2, 2L)
stopifnot(case$wrong == 0, !case$collapsed)

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  replications <- if (length(args)) as.integer(args[1L]) else 16L
  selected <- if (length(args) >= 2L) args[2L] else names(scenarios)
  output <- if (length(args) >= 3L) args[3L] else "benchmarks/stress-results.csv"
  offset <- if (length(args) >= 4L) as.integer(args[4L]) else 0L
  if (!all(is.finite(c(replications, offset))) || replications < 1L ||
      offset < 0L) stop("Invalid replication range")
  rows <- list()
  for (name in selected) {
    if (!name %in% names(scenarios)) stop("Unknown scenario ", name)
    s <- scenarios[[name]]
    for (i in offset + seq_len(replications)) {
      data <- simulate(s, 1000L + i)
      truth <- data$state[-seq_len(data$order)]
      oracle <- list(wrong = NA_real_)
      if (is.null(s$innovation)) {
        oracle <- classify(env$zeta_YT_function(s$beta, data$P, data$y,
          s$k, data$order)[, -seq_len(data$order), drop = FALSE], truth, s$k)
      }
      for (package in c("MSARM", "MSwM")) {
        fit <- fit_one(package, data, s, 2000L + i, return_path = TRUE)
        path <- fit$path
        fit$path <- NULL
        failed <- nzchar(fit$error)
        measured <- if (failed) classify(matrix(1, s$k, length(truth)),
                                         truth, s$k) else classify(path, truth, s$k)
        if (failed) {
          keep <- c("true_states", "true_counts")
          measured[setdiff(names(measured), keep)] <- NA
        }
        rows[[length(rows) + 1L]] <- as.data.frame(c(
          list(scenario = name, replication = i, n_scored = length(truth),
               oracle_wrong = oracle$wrong, no_fit = failed,
               any_error_failure = failed || measured$wrong > 0,
               failure = failed || measured$wrong > 0.2 ||
                 measured$zero_recall > 0L || measured$collapsed),
          fit, measured))
      }
      write.csv(do.call(rbind, rows), output, row.names = FALSE)
      message("Completed ", name, " replication ", i)
    }
  }
}
