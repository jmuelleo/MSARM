# Regime recovery stress test: MSARM 1.0.1 and MSwM 1.5

The test repeats 96 series from the original benchmark and adds 48 series with asymmetric transitions, brief regimes, rapid switching, weak separation, near-unit-root dynamics, and a rare third regime. Each package receives the true AR order, regime count and switching pattern. The fitting budgets match the original benchmark, including five starts for each estimator. The best fit is selected by the same conditional Gaussian log-likelihood.

For each estimate, the same smoother decodes the most probable regime at every post-lag observation. Labels are aligned by the permutation with the fewest errors, with mean recall breaking ties. Wrong time is the percentage of scored observations assigned to the wrong regime. Balanced error averages the error rates of the true regimes that appeared, giving each equal weight. The oracle uses the generating parameters with the same decoder.

A **failure** means no finite fit, more than 20% wrong time, or zero recall for any true regime that appeared. Predicting one regime throughout a multi-regime series therefore always fails, even below 20% wrong time. The strict criterion also counts every fit with any wrongly decoded time point. Classification percentages use successful fits.

## Overall

| Package | Fits | No fit | Single-regime collapse | Any true regime with zero recall | Failure | Strict failure |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| MSARM | 144 | 0 | 0 | 4 | 24 (16.7%) | 143 (99.3%) |
| MSwM | 144 | 4 | 1 | 16 | 75 (52.1%) | 144 (100.0%) |

| Package | All time pooled | Median | 75th pct. | 90th pct. | 95th pct. | Worst | Median balanced error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MSARM | 6.2% | 2.8% | 10.5% | 30.5% | 34.2% | 48.1% | 5.0% |
| MSwM | 17.3% | 11.2% | 38.3% | 45.9% | 48.6% | 55.2% | 25.5% |

| Package | 0% | (0, 5%] | (5%, 10%] | (10%, 20%] | (20%, 30%] | (30%, 50%] | >50% | No fit |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MSARM | 1 | 86 | 16 | 18 | 8 | 15 | 0 | 0 |
| MSwM | 0 | 45 | 19 | 15 | 7 | 52 | 2 | 4 |

## Results by design

| Design | Series | Failures MSARM / MSwM | Collapse MSARM / MSwM | Median wrong time MSARM / MSwM | Oracle median |
| --- | ---: | ---: | ---: | ---: | ---: |
| original | 96 | 14 / 45 | 0 / 0 | 2.5% / 13.0% | 2.4% |
| additional | 48 | 10 / 30 | 0 / 1 | 6.6% / 10.6% | 5.0% |
| intercept | 16 | 0 / 16 | 0 / 0 | 1.6% / 40.2% | 1.4% |
| full | 16 | 0 / 0 | 0 / 0 | 2.1% / 2.0% | 2.0% |
| variance | 16 | 0 / 0 | 0 / 0 | 9.5% / 9.6% | 9.2% |
| mixed_ar2 | 16 | 0 / 6 | 0 / 0 | 2.2% / 2.4% | 1.9% |
| weak_short | 16 | 14 / 15 | 0 / 0 | 34.2% / 37.3% | 21.5% |
| three_ar2 | 16 | 0 / 8 | 0 / 0 | 1.8% / 12.1% | 1.8% |
| rare | 8 | 0 / 7 | 0 / 1 | 0.6% / 4.0% | 0.8% |
| brief | 8 | 1 / 4 | 0 / 0 | 2.4% / 4.0% | 2.0% |
| fast | 8 | 0 / 0 | 0 / 0 | 10.0% / 10.0% | 10.4% |
| weak_long | 8 | 4 / 7 | 0 / 0 | 19.3% / 42.6% | 11.2% |
| near_unit | 8 | 4 / 6 | 0 / 0 | 22.9% / 41.2% | 18.1% |
| rare_three | 8 | 1 / 6 | 0 / 0 | 2.2% / 15.3% | 1.6% |

## Fits that predict a single regime

| Design | Replication | Package | True counts | Predicted counts | Wrong time |
| --- | ---: | --- | ---: | ---: | ---: |
| rare | 2 | MSwM | 242:7 | 249:0 | 2.8% |

The [per-fit results](stress-results.csv) include occupancies and confusion matrices (true regimes in rows, predicted regimes in columns, serialized row by row). [Summary statistics](stress-summary.csv) include results by design. A single true regime in a simulated sample does not count as collapse. Absent true regimes are omitted from per-regime averages.

## Numerical failures and missed regimes

MSwM returned no fit on four rare-third-regime series (replications 2, 3, 5 and 7). Replaying each of the five starts reproduced a singular system in LAPACK's dgesv in every case.

MSARM had zero recall for at least one true regime on four series. In rare_three replication 1, the rare regime occurred once in 248 observations. MSARM assigned seven points to that regime after label alignment, all from another regime, and missed the only true occurrence. The overall error was 4.4%, and the fit counts as a failure.

## Reproduction

From the repository root, use R with MSwM 1.5 installed. The last argument is the replication offset.

```sh
mkdir -p stress-batches
for s in intercept full variance mixed_ar2 weak_short three_ar2; do
  for offset in 0 8; do
    Rscript benchmarks/stress-regimes.R 8 "$s" "stress-batches/$s-$offset.csv" "$offset"
  done
done
for s in rare brief fast weak_long near_unit rare_three; do
  Rscript benchmarks/stress-regimes.R 8 "$s" "stress-batches/$s-0.csv" 0
done
Rscript benchmarks/summarize-stress.R stress-batches
Rscript benchmarks/render-stress.R > benchmarks/REGIME-STRESS.md
```
