# MSARM 1.0.1 versus MSwM 1.5

Both estimators fit the same 96 Gaussian Markov-switching series, with 16
replications of each design. The true order, number of regimes, and switching
pattern were supplied. MSARM used five starts and 250 EM steps per start.
MSwM used five independent fits with a 250-iteration cap and its own
initialization and convergence rules.

Both estimates were scored with the same independent conditional Gaussian
likelihood and a uniform initial regime distribution. Differences of at
most 0.5 log-likelihood units count as ties. Classification error uses the
same smoother for both estimates and the best permutation of the regime
labels against the simulated path.

| Design | MSARM wins | MSwM wins | Ties | Median regime error MSARM / MSwM |
| --- | ---: | ---: | ---: | ---: |
| Intercept switching only | 16 | 0 | 0 | 1.6% / 40.2% |
| All parameters switch | 0 | 0 | 16 | 2.1% / 2.0% |
| Variance switching only | 12 | 0 | 4 | 9.5% / 9.6% |
| AR(2), mixed switching | 6 | 0 | 10 | 2.2% / 2.4% |
| Short series, weak separation | 7 | 2 | 7 | 34.2% / 37.3% |
| Three regimes, AR(2) | 11 | 0 | 5 | 1.8% / 12.1% |
| **Overall** | **52** | **2** | **42** | **2.5% / 13.0%** |

Both packages returned a finite fit for all 96 series. A score more than
two units below the likelihood evaluated at the generating parameters
occurred in 0 MSARM fits and 30 MSwM fits. The paired median likelihood
difference was +0.80 for MSARM. Weakly separated regimes were difficult
to classify with either estimator. These results describe the tested
designs and fitting budgets.

## Reproduction

Run from the repository root, with R and MSwM 1.5 installed.

```sh
Rscript benchmarks/compare-mswm.R 16
Rscript benchmarks/time-original.R
```

[Per-series results](results.csv) include elapsed times. Simulation batches
ran concurrently, so those timings are descriptive. The isolated comparison
against the original MSARM commit used an AR(2) series of 250 observations,
30 EM steps, and two starts. Over three runs median elapsed time was
1.84 seconds for the original and 0.14 seconds for the revision on Windows
with R 4.5.2.
