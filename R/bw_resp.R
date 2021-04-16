#' Butterworth filter frequency response
#'
#' This function calculates Butterworth filter frequency response (for both causal and acausal)
#' @param freq An array of frequency
#' @param fc The corner frequency for low-pass/high-pass
#' @param nPole The pole parameter for the order, where the order is 2 * nPole. Positive for low-pass; Negative for high-pass
#' @param is_causal A binary to indicator if causal is applied. Apply causal if is_causal = TRUE and apply acausal if is_causal = FALSE
#' @return The frequency response
#' @export
bw_resp <- function(freq, fc, nPole, is_causal){

  n <- abs(nPole)
  # causal
  if (is_causal) {
    hs <- complex(real = 1)  # hs = complex filter response at frequency f
    if (n == 0) return(hs)
    if (is.na(fc)) return(hs)
    # if freq = 0 Hz
    ii <- which(freq %in% c(0))
    freq[ii] <- 10^-6
    as <- complex(imaginary = abs(freq / fc))
    if (nPole < 0) as <- 1 / as  # if high-pass
    hs <- as - exp(complex(imaginary = pi * (0.5 + (((2 * 1.) - 1.) / (2. * n)))))
    if (n > 1)
      for (i in 2:n)
        hs <- hs * (as - exp(complex(imaginary = pi * (0.5 + (((2. * i) -1.) / (2. * n))))))
    hs <- 1 / hs
    return(hs)
  }

  # acausal
  if (!is_causal) {
    nNyq <- length(freq)
    deltaf <- freq[3] - freq[2]
    domega <- deltaf / fc
    omega <- rep(0, nNyq)
    omega <- (seq(1:nNyq) - 1) * domega

    if (nPole > 0) {
      lp <- sqrt(1.0 / (1.0 + omega^(2 * n)))
      return(lp)
    }

    if (nPole < 0) {
      hp <- sqrt(omega^(2 * n) / (1.0 + omega^(2 * n)))
      return(hp)
    }
  }
}
