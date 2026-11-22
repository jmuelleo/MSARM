msarm_count <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

msarm_design <- function(y, K) {
  lagged <- stats::embed(as.numeric(y), K + 1L)
  list(y = lagged[, 1L], X = cbind(1, lagged[, -1L, drop = FALSE]))
}

msarm_update <- function(y, X, probabilities, alpha, switching, floor) {
  N <- nrow(alpha)
  p <- ncol(X)
  changing <- switching[seq_len(p)]
  variance <- alpha[, p + 1L]
  mass <- rowSums(probabilities)
  active <- which(mass > 0)
  beta <- alpha[, seq_len(p), drop = FALSE]
  weighted <- probabilities[active, , drop = FALSE] / variance[active]
  fixed <- X[, !changing, drop = FALSE]
  varying <- X[, changing, drop = FALSE]

  # Eliminate the regime-specific columns by weighted QR before solving for
  # shared coefficients. This avoids a large, mostly zero stacked design.
  if (ncol(fixed)) {
    projected <- lapply(seq_along(active), function(i) {
      weight <- sqrt(weighted[i, ])
      block <- cbind(fixed, y) * weight
      if (ncol(varying)) {
        block <- qr.resid(qr(varying * weight), block)
      }
      block
    })
    blocks <- do.call(rbind, projected)
    shared <- stats::lm.fit(blocks[, seq_len(ncol(fixed)), drop = FALSE],
                            blocks[, ncol(blocks)])$coefficients
    shared[is.na(shared)] <- 0
    beta[, !changing] <- rep(shared, each = N)
    y <- y - drop(fixed %*% shared)
  }

  if (ncol(varying)) {
    for (j in active) {
      fit <- stats::lm.wfit(varying, y, probabilities[j, ])
      coef <- fit$coefficients
      coef[is.na(coef)] <- 0
      beta[j, changing] <- coef
    }
  }
  errors <- sweep(varying %*% t(beta[, changing, drop = FALSE]),
                  1L, y, "-")
  ss <- rowSums(probabilities * t(errors^2))
  if (switching[p + 1L]) {
    variance[active] <- pmax(ss[active] / mass[active], floor)
  } else {
    variance[] <- max(sum(ss) / sum(mass), floor)
  }
  cbind(beta, variance)
}

msarm_fit <- function(Y_T, K, N, Switcher, m, threshold, maxit, starts,
                      criterion, all.plot) {
  K <- msarm_count(K, "K")
  N <- msarm_count(N, "N", 2L)
  m <- msarm_count(m, "m")
  maxit <- msarm_count(maxit, "max_value")
  starts <- msarm_count(starts, "R_value")
  criterion <- match.arg(criterion, c("LV", "RSS", "RCM", "Entropy"))
  if (!is.numeric(Y_T) || !is.null(dim(Y_T)) || any(!is.finite(Y_T)) ||
      length(Y_T) - K < max(K + 2, m + 1)) {
    stop("Y_T must be a finite univariate series with enough observations.",
         call. = FALSE)
  }
  if (!is.logical(Switcher) || length(Switcher) != K + 2L ||
      anyNA(Switcher)) {
    stop("Switcher must contain K + 2 logical values.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0 || threshold > 1 ||
      !is.logical(all.plot) || length(all.plot) != 1L || is.na(all.plot)) {
    stop("Supply a threshold in [0, 1] and a logical all.plot value.",
         call. = FALSE)
  }

  design <- msarm_design(Y_T, K)
  floor <- stats::var(Y_T) * 1e-10
  if (!is.finite(floor) || floor <= 0) {
    stop("Y_T must have positive, finite variance.", call. = FALSE)
  }
  times <- seq.int(K + 1L, length(Y_T))
  fits <- vector("list", starts)
  for (r in seq_len(starts)) {
    alpha <- alpha_0_function(Y_T, N, K)
    P <- P_0_function(Y_T, N, K)
    filter <- msarm_filter(alpha, P, Y_T, N, K, design)
    for (iteration in seq_len(maxit)) {
      smooth <- msarm_smooth(filter, P, K)
      alpha <- msarm_update(design$y, design$X, smooth[, times],
                            alpha, Switcher, floor)
      P <- msarm_transition(P, filter, smooth, K, m)
      filter <- msarm_filter(alpha, P, Y_T, N, K, design)
    }
    smooth <- msarm_smooth(filter, P, K)[, times]
    regime <- max.col(t(smooth), ties.method = "first")
    fitted <- rowSums(design$X * alpha[regime, seq_len(K + 1L), drop = FALSE])
    residuals <- design$y - fitted
    positive <- smooth[smooth > 0]
    fits[[r]] <- list(P = P, alpha = alpha, zeta_t_T = smooth,
                      LV = filter$loglik, RSS = sum(residuals^2),
                      RCM = 100 * N / (N - 1) * mean(colSums(smooth * (1 - smooth))),
                      Entropy = -sum(positive * log(positive)) /
                        (ncol(smooth) * log(N)),
                      Y_T = Y_T, Yhat_insample = fitted, residuals = residuals)
    if (all.plot) {
      graphics::matplot(t(smooth), type = "l", lty = 1,
                        xlab = "Observation", ylab = "Regime probability")
      graphics::abline(h = threshold, lty = 2)
      graphics::matplot(cbind(design$y, fitted), type = "l", lty = 1,
                        xlab = "Observation", ylab = "Observed and fitted")
    }
  }
  scores <- vapply(fits, `[[`, numeric(1), criterion)
  best <- if (criterion == "LV") which.max(scores) else which.min(scores)
  fits[[best]]
}

msarm_update_from_series <- function(Y_T, alpha, P, Switcher, N, K) {
  design <- msarm_design(Y_T, K)
  smooth <- zeta_YT_function(alpha, P, Y_T, N, K)
  msarm_update(design$y, design$X, smooth[, -seq_len(K), drop = FALSE],
               alpha, Switcher, stats::var(Y_T) * 1e-10)
}
