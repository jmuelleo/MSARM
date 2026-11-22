# MSARM 1.0.1

- Fit all switching patterns with one Gaussian EM implementation. The
  likelihood and the coefficient update use the full post-lag sample.
  A uniform initial regime distribution makes likelihood evaluations
  deterministic. Densities are evaluated in log space.
- Use weighted QR for shared and regime-specific coefficients. The variance
  update uses a lower bound equal to 1e-10 times the sample variance.
- Compute conditional multi-step forecasts and draw entire residual
  bootstrap paths. Forecast plots include the forecast range and retain the
  time index of `ts` inputs by default.
- Keep the existing function names and result fields. The transition matrix
  remains column-stochastic. The classification measure now uses the
  0-to-100 generalization for multiple regimes. Entropy is normalized over
  all regimes. Fitting is silent unless plots are requested. Presample
  entries of probability and density matrices are now marked with `NA`.
