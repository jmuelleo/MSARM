# MSARM

MSARM fits Gaussian Markov-switching autoregressions. The intercept, each
lag coefficient, and the innovation variance can switch independently.
The original thesis describes the model and its derivation.
[Theory and implementation](https://github.com/jmuelleo/MSARM/blob/master/MSARM_Theory_and_Implementation.pdf)

```r
set.seed(12)
y <- arima.sim(list(ar = 0.5), n = 150)
fit <- MSARM::MSARM.fit(y, K = 1, N = 2,
                        Switcher = c(TRUE, FALSE, TRUE))
forecast <- MSARM::MSARM.predict(fit, n.ahead = 8, boot = TRUE, L = 1000)
MSARM::MSARM.plot(forecast, conf = TRUE)
```

`fit$alpha` has one row per regime and columns for the intercept, AR
coefficients and innovation variance. `fit$P[j, i]` is the probability
of moving from regime `i` to regime `j`. Columns of `fit$zeta_t_T`
are smoothed probabilities for observations `K + 1` through
`length(y)`. Regime labels are arbitrary.

The fitted regression assumes Gaussian innovations conditional on the regime.
With switching variances the likelihood can favor near-degenerate regimes.
The estimator applies a variance floor and uses multiple random starts.
Use `set.seed()` to reproduce a fit.

[Reproducible comparison with MSwM](https://github.com/jmuelleo/MSARM/blob/msarm-1.0.1/benchmarks/README.md)
