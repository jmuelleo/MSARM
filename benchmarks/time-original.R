# Run from the repository root with Rscript benchmarks/time-original.R.
revision <- "f1aa60cfa676b06ca7628a31d0ce9e1d1629bc5c"
original <- system2("git", c("ls-tree", "-r", "--name-only", revision, "R"),
                    stdout = TRUE)
if (!length(original) || !is.null(attr(original, "status"))) {
  stop("The original MSARM commit is required.")
}
old <- new.env()
new <- new.env()
for (file in original) {
  text <- system2("git", c("show", paste0(revision, ":", file)), stdout = TRUE)
  eval(parse(text = text), envir = old)
}
invisible(lapply(list.files("R", full.names = TRUE), sys.source, envir = new))
set.seed(18)
y <- as.numeric(arima.sim(list(ar = c(0.5, -0.15)), n = 250))
timings <- function(env) {
  replicate(3, system.time({
    set.seed(123)
    invisible(capture.output(env$MSARM.fit(
      y, 2, 2, c(TRUE, FALSE, TRUE, TRUE), max_value = 30, R_value = 2)))
  })[["elapsed"]])
}
pdf(NULL)
results <- rbind(original = timings(old), revised = timings(new))
invisible(dev.off())
print(results)
print(apply(results, 1L, median))
