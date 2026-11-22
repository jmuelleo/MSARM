# Run from the repository root with the directory of 26 new batch CSVs.
batch_dir <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(batch_dir)) stop("Supply the directory of the 26 new batch CSVs.")
source("benchmarks/stress-regimes.R")
first <- read.csv("benchmarks/stress-results.csv")
existing <- unique(first$scenario)
new <- setdiff(names(scenarios), existing)
files <- file.path(batch_dir, paste0(rep(new, each = 2L), "-", c(0, 6), ".csv"))
if (!all(file.exists(files))) stop("Missing simulation batches.")
batches <- lapply(files, read.csv)
stopifnot(length(new) == 13L, all(vapply(batches, nrow, 0L) == 12L))
fits <- rbind(first, do.call(rbind, batches))
fits <- fits[order(match(fits$scenario, names(scenarios)), fits$replication,
                   fits$package), ]
key <- function(x) paste(x$scenario, x$replication, sep = "/")
stopifnot(nrow(fits) == 600L, length(unique(key(fits))) == 300L,
          all(table(key(fits), fits$package) == 1L),
          all(fits$no_fit == is.na(fits$wrong)),
          isTRUE(all.equal(fits$wrong, fits$misclassified)))
pairs <- fits[fits$package == "MSARM", c("scenario", "replication")]
series <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  row <- pairs[i, ]
  s <- scenarios[[row$scenario]]
  data <- simulate(s, 1000L + row$replication)
  truth <- data$state[-seq_len(data$order)]
  counts <- tabulate(truth, s$k)
  generating_ll <- NA_real_
  if (is.null(s$innovation)) {
    generating_ll <- score(data$y, data$order, data$P, s$beta)
  }
  data.frame(row, k = s$k, order = data$order,
    n_scored = length(truth), true_states = sum(counts > 0),
    true_counts = paste(counts, collapse = ":"),
    min_prevalence = min(counts[counts > 0]) / length(truth),
    switches = sum(diff(truth) != 0),
    generating_ll = generating_ll)
}))
paired <- match(key(fits), key(series))
stopifnot(identical(fits$true_counts, series$true_counts[paired]),
          all(fits$n_scored == series$n_scored[paired]))
fits$min_prevalence <- series$min_prevalence[paired]
fits$true_switches <- series$switches[paired]
fits$k <- series$k[paired]
fits$order <- series$order[paired]
fits$ll_gap <- fits$loglik - series$generating_ll[paired]
fits$oracle_excess <- fits$wrong - fits$oracle_wrong
fits$spurious_split <- with(fits, true_states == 1 & predicted_states > 1)
fits$failure_reason <- with(fits, ifelse(no_fit, "no fit",
  ifelse(collapsed, "single-regime collapse", ifelse(zero_recall > 0,
    "other zero recall", ifelse(wrong > .2, "over 20% wrong", "passed")))))
stopifnot(all(fits$failure == (fits$failure_reason != "passed")),
          all(vapply(which(!fits$no_fit), function(i) {
            x <- fits[i, ]
            read_counts <- function(s) as.integer(strsplit(s, ":", fixed = TRUE)[[1L]])
            m <- matrix(read_counts(x$confusion), x$k, byrow = TRUE)
            all(rowSums(m) == read_counts(x$true_counts)) &&
              all(colSums(m) == read_counts(x$predicted_counts)) &&
              abs(1 - sum(diag(m)) / sum(m) - x$wrong) < 1e-12
          }, TRUE)))

stats <- function(x) {
  ok <- x[!x$no_fit, ]
  w <- ok$wrong
  n <- nrow(x)
  q <- function(p) if (length(w)) quantile(w, p, names = FALSE) else NA_real_
  count <- function(z) sum(z, na.rm = TRUE)
  fails_at <- function(cutoff) count(x$no_fit | x$zero_recall > 0L | x$wrong > cutoff)
  c(n = n, failure = count(x$failure), no_fit = count(x$no_fit),
    fail_05 = fails_at(.05), fail_10 = fails_at(.10), fail_30 = fails_at(.30),
    fail_50 = fails_at(.50),
    zero_recall = count(ok$zero_recall > 0L), collapse = count(ok$collapsed),
    over20 = count(w > .2),
    only_one_true = count(x$true_states == 1L),
    spurious_split = count(ok$spurious_split),
    strict_failure = count(x$any_error_failure),
    pooled_wrong = weighted.mean(w, ok$n_scored),
    mean_wrong = mean(w), median = q(.5), p75 = q(.75), p90 = q(.9),
    p95 = q(.95), worst = q(1),
    median_balanced = median(ok$balanced_wrong),
    median_min_recall = median(ok$min_recall),
    median_oracle = median(x$oracle_wrong, na.rm = TRUE),
    median_oracle_excess = median(ok$oracle_excess, na.rm = TRUE),
    above_oracle_10pp = count(ok$oracle_excess > .10),
    likelihood_below_truth_2 = count(ok$ll_gap < -2),
    exact = count(w == 0), up_to5 = count(w > 0 & w <= .05),
    up_to10 = count(w > .05 & w <= .10), up_to20 = count(w > .10 & w <= .20),
    up_to30 = count(w > .20 & w <= .30), up_to50 = count(w > .30 & w <= .50),
    above50 = count(w > .50),
    elapsed_median = median(ok$seconds))
}
groups <- c(list(all = rep(TRUE, nrow(fits)),
                 previous = fits$scenario %in% existing,
                 new = fits$scenario %in% new,
                 observed_rare = fits$true_states >= 2 & fits$min_prevalence < .1,
                 one_observed_regime = fits$true_states == 1,
                 three_regimes = fits$k == 3,
                 short_n_under_100 = fits$n_scored < 100,
                 misspecified = fits$scenario %in% c("heavy_tail", "outliers")),
            setNames(lapply(names(scenarios), function(s) fits$scenario == s),
                     names(scenarios)))
summary <- do.call(rbind, lapply(names(groups), function(g) {
  do.call(rbind, lapply(c("MSARM", "MSwM"), function(p) {
    x <- fits[groups[[g]] & fits$package == p, ]
    data.frame(group = g, package = p, as.list(stats(x)))
  }))
}))
stopifnot(all(with(summary,
  exact + up_to5 + up_to10 + up_to20 + up_to30 + up_to50 + above50 + no_fit == n)))
write.csv(fits, "benchmarks/stress300-results.csv", row.names = FALSE)
write.csv(summary, "benchmarks/stress300-summary.csv", row.names = FALSE)
write.csv(series, "benchmarks/stress300-series.csv", row.names = FALSE)
print(summary[, c("group", "package", "n", "failure", "no_fit", "zero_recall",
                  "collapse", "median", "median_oracle")], row.names = FALSE)
