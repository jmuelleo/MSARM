# Regime recovery on 300 paired simulated series

This study extends the [144-series stress test](REGIME-STRESS.md) with 156 new series, 12 per new design. There are 25 fixed designs in total. The comparisons use MSARM 1.0.1 and MSwM 1.5, with the same series, true AR order, number of regimes and switching mask supplied to each package. MSARM uses five starts and 250 EM iterations per start. MSwM is called five times with a 250-iteration cap, retaining its own initialization and convergence criteria. Candidates are scored using the same conditional Gaussian likelihood. For both packages, classification uses the same smoother, posterior mode and permutation of regime labels that minimizes the number of wrong assignments.

A failure means no finite fit, over 20% wrong time, or zero recall for any observed regime. Collapse means one regime predicted at every timestamp on a series with multiple observed regimes. It is a subset of zero recall, and always fails. The first AR lags are excluded. Pooled wrong time weights each successful fit by the number of scored observations, while quantiles weight series equally. The oracle applies the shared decoder to the true Gaussian model. It is a reference for classification difficulty, and may be outperformed on individual finite samples. It is omitted for the Student-t and contaminated Gaussian designs.

## Overall results

| Outcome | MSARM | MSwM |
| --- | ---: | ---: |
| Failure | 78/300 (26.0%) | 189/300 (63.0%) |
| No finite fit | 0 | 4 |
| At least one observed regime with zero recall | 10 | 49 |
| Single-regime collapse | 0 | 5 |
| Wrong time, pooled | 9.5% | 22.2% |
| Wrong time, median per series | 3.1% | 24.4% |
| 75th percentile | 18.7% | 41.1% |
| 90th percentile | 34.2% | 48.3% |
| 95th percentile | 44.2% | 49.4% |
| Worst series | 49.7% | 61.0% |
| Median equal-weight error across observed regimes | 5.5% | 37.1% |
| Any wrong time point or no fit | 241/300 | 283/300 |

Among the 296 pairs with two finite fits, MSARM classified more timestamps correctly in 176 cases and MSwM in 45, with 75 ties. Only MSARM failed on 5 series, only MSwM on 116, both on 73, and neither on 106. On the new 156 series alone, the failure counts were 54 and 114.

| Per-series wrong time | MSARM | MSwM |
| --- | ---: | ---: |
| Exactly 0% | 59 | 17 |
| Greater than 0%, up to 5% | 106 | 54 |
| Greater than 5%, up to 10% | 28 | 34 |
| Greater than 10%, up to 20% | 33 | 25 |
| Greater than 20%, up to 30% | 33 | 37 |
| Greater than 30%, up to 50% | 41 | 123 |
| Greater than 50% | 0 | 6 |
| No fit | 0 | 4 |

The advantage persists across the tested error cutoffs, while always counting zero recall and no fit.

| Wrong-time cutoff | 5% | 10% | 20% | 30% | 50% |
| --- | ---: | ---: | ---: | ---: | ---: |
| MSARM failures | 136 | 109 | 78 | 47 | 10 |
| MSwM failures | 240 | 213 | 189 | 161 | 57 |

## Where results diverge

- **Well-separated intercepts.** In the new clear two-regime, short clear, and clear three-regime scenarios, MSARM had 0 failures in each of the three 12-series groups. MSwM had 10, 7 and 9. For `clear/1`, the true intercepts are -3 and +3 with AR coefficient 0.20. All five MSwM starts reached likelihoods near -293.75, with nearly equal intercepts and an AR coefficient around 0.84, essentially equal to the series' single-regime AR estimate of 0.8431. The best MSwM solution assigned all 159 timestamps to one regime. MSARM recovered the partition exactly at log-likelihood -196.11. Running 250 MSARM updates from the best MSwM parameters barely moved the likelihood, to -293.75, so changing the update rule alone did not recover the high-likelihood solution from that starting point.
- **Full switching and clear dynamics.** Both packages passed all 16 original full-switching series and all 8 fast-switching series with changes in multiple parameters. With AR coefficients alone switching, both failed on 2 of 12, with identical 8.4% median wrong time. The difference is therefore concentrated in particular designs.
- **Weak signal.** On the strongly overlapping design the generating-parameter classifier already missed a median 30.5% of timestamps, and MSARM and MSwM failed on 12 and 11 of 12. With fast variance-only switching both failed on all 12 and the oracle median was 22.0%. For `quick_var/1`, the generating decoder missed 22.0%, MSARM missed 23.9%, and the best original MSwM fit missed 27.0%. Re-estimating from the MSwM solution with MSARM converged to the same 23.9% solution as the independent MSARM fit. More iterations alone cannot remove classification uncertainty caused by overlapping emissions.
- **Rare regimes.** Of 48 series with a visited regime below 10% prevalence and multiple visited regimes, MSARM failed on 11 and MSwM on 32. In the rare-and-weak design alone, both failed 10 of 12. MSwM achieved lower median wrong time (10.7% vs 26.4%) while having zero recall for a regime in more series (8 vs 4). The pooled or median time error can favor a model that mostly predicts the majority regime. On three near-absorbing series only one regime actually appeared. Both fitted a second predicted state in all three, and both crossed the 20% error threshold on two.
- **Misspecification.** The standardized Student-t(3) design yielded 0/12 failures for MSARM and 12/12 for MSwM. With Gaussian innovations inflated eightfold at 2.5% of draws, failures were 7/12 and 12/12. These experiments retain Gaussian estimation in both packages. The two designs do not identify the isolated causal effect of tail thickness.
- **Numerical versus statistical failure.** All four MSwM no-fits occurred in the rare-three design. The 144-series follow-up reproduced LAPACK singular systems in every start for all four series. Finite fits can also be poor solutions: on the Gaussian series, 92 MSwM fits and zero MSARM fits scored more than two log-likelihood units below the generating parameter vector. A high score can coexist with poor regime recovery or overfitting. In `rare_three/2`, for example, MSARM scored 42.5 units above the generating vector but drove one variance to about 2.50e-10 for a regime assigned only three points. This is a near-degenerate estimate even though it passed the classification criterion.

Replaying all 300 MSARM fits reproduced their recorded likelihoods and found three variance-floor hits. The affected regimes had effective sample sizes of approximately 2, 3 and 3. One case already failed the classification rule. Counting the two additional fits as failures for numerical degeneracy would put MSARM at 80/300 (26.7%). This supplemental check applies to the variance floor in MSARM.

## By design

Counts are MSARM / MSwM. The first eight rows are aggregates, which can overlap.

| Design or stratum | Series | Failures | No fit | Zero recall | Median wrong time | Oracle median |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| all | 300 | 78 / 189 | 0 / 4 | 10 / 49 | 3.1% / 24.4% | 2.6% |
| previous | 144 | 24 / 75 | 0 / 4 | 4 / 16 | 2.8% / 11.2% | 2.5% |
| new | 156 | 54 / 114 | 0 / 0 | 6 / 33 | 5.0% / 28.4% | 3.1% |
| observed_rare | 48 | 11 / 32 | 0 / 4 | 5 / 22 | 1.2% / 4.5% | 1.0% |
| one_observed_regime | 3 | 2 / 2 | 0 / 0 | 0 / 0 | 23.0% / 23.0% | 0.0% |
| three_regimes | 36 | 1 / 23 | 0 / 4 | 1 / 6 | 1.5% / 23.3% | 1.0% |
| short_n_under_100 | 40 | 18 / 26 | 0 / 0 | 2 / 5 | 11.6% / 34.8% | 9.5% |
| misspecified | 24 | 7 / 24 | 0 / 0 | 1 / 6 | 2.5% / 40.5% | - |
| intercept | 16 | 0 / 16 | 0 / 0 | 0 / 4 | 1.6% / 40.2% | 1.4% |
| full | 16 | 0 / 0 | 0 / 0 | 0 / 0 | 2.1% / 2.0% | 2.0% |
| variance | 16 | 0 / 0 | 0 / 0 | 0 / 0 | 9.5% / 9.6% | 9.2% |
| mixed_ar2 | 16 | 0 / 6 | 0 / 0 | 0 / 0 | 2.2% / 2.4% | 1.9% |
| weak_short | 16 | 14 / 15 | 0 / 0 | 2 / 1 | 34.2% / 37.3% | 21.5% |
| three_ar2 | 16 | 0 / 8 | 0 / 0 | 0 / 0 | 1.8% / 12.1% | 1.8% |
| rare | 8 | 0 / 7 | 0 / 0 | 0 / 6 | 0.6% / 4.0% | 0.8% |
| brief | 8 | 1 / 4 | 0 / 0 | 0 / 4 | 2.4% / 4.0% | 2.0% |
| fast | 8 | 0 / 0 | 0 / 0 | 0 / 0 | 10.0% / 10.0% | 10.4% |
| weak_long | 8 | 4 / 7 | 0 / 0 | 0 / 0 | 19.3% / 42.6% | 11.2% |
| near_unit | 8 | 4 / 6 | 0 / 0 | 1 / 0 | 22.9% / 41.2% | 18.1% |
| rare_three | 8 | 1 / 6 | 0 / 4 | 1 / 1 | 2.2% / 15.3% | 1.6% |
| clear | 12 | 0 / 10 | 0 / 0 | 0 / 4 | 0.0% / 33.0% | 0.0% |
| rare_clear | 12 | 0 / 4 | 0 / 0 | 0 / 3 | 0.0% / 1.1% | 0.0% |
| rare_weak | 12 | 10 / 10 | 0 / 0 | 4 / 8 | 26.4% / 10.7% | 4.4% |
| overlap | 12 | 12 / 11 | 0 / 0 | 1 / 2 | 43.0% / 46.0% | 30.5% |
| ar_only | 12 | 2 / 2 | 0 / 0 | 0 / 0 | 8.4% / 8.4% | 7.5% |
| short_clear | 12 | 0 / 7 | 0 / 0 | 0 / 2 | 0.0% / 31.5% | 0.0% |
| short_var | 12 | 4 / 4 | 0 / 0 | 0 / 2 | 10.9% / 14.5% | 8.7% |
| quick_var | 12 | 12 / 12 | 0 / 0 | 0 / 0 | 26.1% / 26.4% | 22.0% |
| clear_three | 12 | 0 / 9 | 0 / 0 | 0 / 5 | 0.0% / 48.7% | 0.0% |
| near_absorb | 12 | 2 / 10 | 0 / 0 | 0 / 1 | 0.0% / 30.5% | 0.0% |
| heavy_tail | 12 | 0 / 12 | 0 / 0 | 0 / 5 | 1.1% / 38.3% | - |
| outliers | 12 | 7 / 12 | 0 / 0 | 1 / 1 | 30.2% / 42.2% | - |
| long_weak | 12 | 5 / 11 | 0 / 0 | 0 / 0 | 16.8% / 45.4% | 14.7% |

## Reproduction and files

The original 144 results are read from [stress-results.csv](stress-results.csv). Reproduce those batches with the instructions in [REGIME-STRESS.md](REGIME-STRESS.md). From the repository root, run the following with MSwM 1.5 installed.

```sh
mkdir -p stress300-batches
for s in clear rare_clear rare_weak overlap ar_only short_clear short_var quick_var clear_three near_absorb heavy_tail outliers long_weak; do
  for offset in 0 6; do
    Rscript benchmarks/stress-regimes.R 6 "$s" "stress300-batches/$s-$offset.csv" "$offset"
  done
done
Rscript benchmarks/analyze-stress300.R stress300-batches
Rscript benchmarks/diagnose-stress300.R diagnostics.csv clear/1 intercept/1 mixed_ar2/1 rare_three/2 rare/2 quick_var/1 heavy_tail/1 outliers/1
```

- [Per-fit results and confusion matrices](stress300-results.csv)
- [Summary by scenario and thresholds](stress300-summary.csv)
- [Series characteristics](stress300-series.csv)
- [Individual diagnostic fits](stress300-diagnostics.csv)
- [MSARM variance audit](stress300-variance-audit.csv). Run the audit with batch numbers 1 through 12, using a separate output CSV for each batch.

The empirical rates apply to this fixed mix of simulated designs and the supplied model specification and fitting budgets. They do not measure how often each design occurs in applications. Elapsed times in the raw data record concurrent runs and do not provide an isolated timing comparison.
