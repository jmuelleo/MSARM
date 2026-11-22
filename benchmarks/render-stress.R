# Run after summarize-stress.R from the repository root.
summary <- read.csv("benchmarks/stress-summary.csv")
fits <- read.csv("benchmarks/stress-results.csv")
cat <- function(..., sep = "") base::cat(..., sep = sep)
pct <- function(x) sprintf("%.1f%%", 100 * x)
line <- function(...) paste0("| ", paste(..., sep = " | "), " |")
cat("# Regime recovery stress test: MSARM 1.0.1 and MSwM 1.5\n\n",
    "The test repeats 96 series from the original benchmark and adds 48 series ",
    "with asymmetric transitions, brief regimes, rapid switching, weak separation, ",
    "near-unit-root dynamics, and a rare third regime. Each package receives the ",
    "true AR order, regime count and switching pattern. The fitting budgets match ",
    "the original benchmark, including five starts for each estimator. The best ",
    "fit is selected by the same conditional Gaussian log-likelihood.\n\n",
    "For each estimate, the same smoother decodes the most probable regime at ",
    "every post-lag observation. Labels are aligned by the permutation with the ",
    "fewest errors, with mean recall breaking ties. Wrong time is the percentage ",
    "of scored observations assigned to the wrong regime. Balanced error averages ",
    "the error rates of the true regimes that appeared, giving each equal weight. ",
    "The oracle uses the generating parameters with the same decoder.\n\n",
    "A **failure** means no finite fit, more than 20% wrong time, or zero recall ",
    "for any true regime that appeared. Predicting one regime throughout a ",
    "multi-regime series therefore always fails, even below 20% wrong time. ",
    "The strict criterion also counts every fit with any wrongly decoded time point. ",
    "Classification percentages use successful fits.\n\n",
    "## Overall\n\n", sep = "")
cat(line("Package", "Fits", "No fit", "Single-regime collapse",
         "Any true regime with zero recall", "Failure", "Strict failure"), "\n")
cat(line("---", "---:", "---:", "---:", "---:", "---:", "---:"), "\n")
for (i in which(summary$group == "all")) {
  x <- summary[i, ]
  cat(line(x$package, x$n, x$no_fit, x$collapsed, x$zero_recall,
           paste0(x$failure, " (", pct(x$failure / x$n), ")"),
           paste0(x$any_error_failure, " (", pct(x$any_error_failure / x$n), ")")),
      "\n")
}
cat("\n")
cat(line("Package", "All time pooled", "Median", "75th pct.", "90th pct.",
         "95th pct.", "Worst", "Median balanced error"), "\n")
cat(line("---", "---:", "---:", "---:", "---:", "---:", "---:", "---:"), "\n")
for (i in which(summary$group == "all")) {
  x <- summary[i, ]
  cat(line(x$package, pct(x$pool_wrong), pct(x$median), pct(x$q75), pct(x$q90),
           pct(x$q95), pct(x$max_wrong), pct(x$median_balanced_wrong)), "\n")
}
cat("\n")
cat(line("Package", "0%", "(0, 5%]", "(5%, 10%]", "(10%, 20%]",
         "(20%, 30%]", "(30%, 50%]", ">50%", "No fit"), "\n")
cat(line("---", "---:", "---:", "---:", "---:", "---:", "---:", "---:", "---:"), "\n")
for (i in which(summary$group == "all")) {
  x <- summary[i, ]
  cat(line(x$package, x$correct, x$wrong_0_to_5, x$wrong_5_to_10,
           x$wrong_10_to_20, x$wrong_20_to_30, x$wrong_30_to_50,
           x$wrong_over50, x$no_fit), "\n")
}
cat("\n## Results by design\n\n")
cat(line("Design", "Series", "Failures MSARM / MSwM", "Collapse MSARM / MSwM",
         "Median wrong time MSARM / MSwM", "Oracle median"), "\n")
cat(line("---", "---:", "---:", "---:", "---:", "---:"), "\n")
for (group in unique(summary$group)[-1L]) {
  x <- summary[summary$group == group, ]
  cat(line(group, x$n[1L], paste(x$failure, collapse = " / "),
           paste(x$collapsed, collapse = " / "), paste(pct(x$median), collapse = " / "),
           pct(x$median_oracle_wrong[1L])), "\n")
}
events <- fits[!is.na(fits$collapsed) & fits$collapsed, ]
cat("\n## Fits that predict a single regime\n\n")
if (nrow(events)) {
  cat(line("Design", "Replication", "Package", "True counts", "Predicted counts",
           "Wrong time"), "\n")
  cat(line("---", "---:", "---", "---:", "---:", "---:"), "\n")
  for (i in seq_len(nrow(events))) {
    x <- events[i, ]
    cat(line(x$scenario, x$replication, x$package, x$true_counts,
             x$predicted_counts, pct(x$wrong)), "\n")
  }
} else {
  cat("No fitted series met the single-regime collapse condition.\n")
}
cat("\nThe [per-fit results](stress-results.csv) include occupancies and confusion ",
    "matrices (true regimes in rows, predicted regimes in columns, serialized ",
    "row by row). [Summary statistics](stress-summary.csv) include results by ",
    "design. A single true regime in a simulated sample does not count as ",
    "collapse. Absent true regimes are omitted from per-regime averages.\n", sep = "")
cat("\n## Numerical failures and missed regimes\n\n",
    "MSwM returned no fit on four rare-third-regime series (replications 2, 3, ",
    "5 and 7). Replaying each of the five starts reproduced a singular ",
    "system in LAPACK's dgesv in every case.\n\n",
    "MSARM had zero recall for at least one true regime on four series. ",
    "In rare_three replication 1, the rare regime occurred once in 248 ",
    "observations. MSARM assigned seven points to that regime after label ",
    "alignment, all from another regime, and missed the only true occurrence. ",
    "The overall error was 4.4%, and the fit counts as a failure.\n\n",
    "## Reproduction\n\n",
    "From the repository root, use R with MSwM 1.5 installed. ",
    "The last argument is the replication offset.\n\n",
    "```sh\n",
    "mkdir -p stress-batches\n",
    "for s in intercept full variance mixed_ar2 weak_short three_ar2; do\n",
    "  for offset in 0 8; do\n",
    '    Rscript benchmarks/stress-regimes.R 8 "$s" "stress-batches/$s-$offset.csv" "$offset"\n',
    "  done\n",
    "done\n",
    "for s in rare brief fast weak_long near_unit rare_three; do\n",
    '  Rscript benchmarks/stress-regimes.R 8 "$s" "stress-batches/$s-0.csv" 0\n',
    "done\n",
    "Rscript benchmarks/summarize-stress.R stress-batches\n",
    "Rscript benchmarks/render-stress.R > benchmarks/REGIME-STRESS.md\n",
    "```\n", sep = "")
