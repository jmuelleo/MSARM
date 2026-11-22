# Run from the repository root with the directory of per-batch CSVs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the directory of per-batch CSVs.")
original <- c("intercept", "full", "variance", "mixed_ar2", "weak_short", "three_ar2")
extra <- c("rare", "brief", "fast", "weak_long", "near_unit", "rare_three")
designs <- c(original, extra)
expected <- c(outer(original, c(0, 8), paste, sep = "-"),
              paste0(extra, "-0"))
files <- file.path(args[1L], paste0(expected, ".csv"))
if (!all(file.exists(files))) stop("Missing simulation batches.")
batches <- lapply(files, read.csv, stringsAsFactors = FALSE)
stopifnot(all(lengths(lapply(batches, rownames)) == 16L))
fits <- do.call(rbind, batches)
fits <- fits[order(match(fits$scenario, designs), fits$replication,
                   fits$package), ]
stopifnot(!anyDuplicated(fits[c("scenario", "replication", "package")]),
          nrow(fits) == 288L,
          all(fits$no_fit == is.na(fits$wrong)),
          isTRUE(all.equal(fits$wrong, fits$misclassified)),
          all(fits$failure[!is.na(fits$collapsed) & fits$collapsed]))
check_counts <- function(i) {
  split_counts <- function(field) as.integer(strsplit(field, ":", fixed = TRUE)[[1L]])
  row <- fits[i, ]
  true <- split_counts(row$true_counts)
  predicted <- split_counts(row$predicted_counts)
  confusion <- matrix(split_counts(row$confusion), length(true), byrow = TRUE)
  identical(as.numeric(true), rowSums(confusion)) &&
    identical(as.numeric(predicted), colSums(confusion)) &&
    sum(true) == row$n_scored &&
    abs(1 - sum(diag(confusion)) / sum(true) - row$wrong) < 1e-12
}
stopifnot(all(vapply(which(!fits$no_fit), check_counts, TRUE)))
write.csv(fits, "benchmarks/stress-results.csv", row.names = FALSE)

summarize <- function(x) {
  ok <- x[!x$no_fit, ]
  wrong <- ok$wrong
  count <- function(x) sum(x, na.rm = TRUE)
  missed <- ok$zero_recall > 0L
  c(n = nrow(x), no_fit = count(x$no_fit),
    collapsed = count(ok$collapsed),
    never_predicted = count(ok$never_predicted > 0L),
    zero_recall = count(missed),
    zero_recall_noncollapse = count(missed & !ok$collapsed),
    over20_without_missed = count(wrong > .2 & !missed & !ok$collapsed),
    failure = count(x$failure),
    any_error_failure = count(x$any_error_failure),
    correct = count(wrong == 0),
    wrong_0_to_5 = count(wrong > 0 & wrong <= .05),
    wrong_5_to_10 = count(wrong > .05 & wrong <= .10),
    wrong_10_to_20 = count(wrong > .10 & wrong <= .20),
    wrong_20_to_30 = count(wrong > .20 & wrong <= .30),
    wrong_30_to_50 = count(wrong > .30 & wrong <= .50),
    wrong_over50 = count(wrong > .50),
    pool_wrong = weighted.mean(wrong, ok$n_scored),
    mean_wrong = mean(wrong),
    q10 = quantile(wrong, .10, names = FALSE),
    median = median(wrong),
    q75 = quantile(wrong, .75, names = FALSE),
    q90 = quantile(wrong, .90, names = FALSE),
    q95 = quantile(wrong, .95, names = FALSE),
    max_wrong = max(wrong),
    median_balanced_wrong = median(ok$balanced_wrong),
    median_min_recall = median(ok$min_recall),
    median_oracle_wrong = median(x$oracle_wrong),
    one_true_state = count(x$true_states == 1L))
}
groups <- c(list(all = designs, original = original, additional = extra),
            as.list(setNames(designs, designs)))
summary <- do.call(rbind, lapply(names(groups), function(label) {
  cases <- fits[fits$scenario %in% groups[[label]], ]
  do.call(rbind, lapply(c("MSARM", "MSwM"), function(package) {
    stats <- summarize(cases[cases$package == package, ])
    data.frame(group = label, package = package, as.list(stats))
  }))
}))
write.csv(summary, "benchmarks/stress-summary.csv", row.names = FALSE)
print(summary[summary$group == "all", ], row.names = FALSE)
