#'@title hmm_genotyping
#'@description Infer haplotype inheritance based on Hidden Markov Model
hmm_genotyping <- function(vcfdat, genetic_map) {
  initial_states <- c(0.5, 0.5)
  symbols <- seq_len(nrow(vcfdat))
  observation <- seq_len(nrow(vcfdat))
  emissionProbs <- t(vcfdat[, c("P0", "P1")])
  transProbs <- getTransProbs(genetic_map, vcfdat)

  # construct HMM model
  hmm <-
    initHMM(c("1", "2"),
            symbols,
            initial_states,
            transProbs,
            emissionProbs)
  vcfdat[, hmm_fetal_hap := as.numeric(viterbi(hmm, observation))]

  vcfdat$P0re <-
    c(0, sapply(seq_len(nrow(vcfdat) - 1), function(n) {
      prev_state = vcfdat$hmm_fetal_hap[n]
      return(transProbs[prev_state, 1, n])
    }))
  vcfdat$P1re <-
    c(0, sapply(seq_len(nrow(vcfdat) - 1), function(n) {
      prev_state = vcfdat$hmm_fetal_hap[n]
      return(transProbs[prev_state, 2, n])
    }))
  
  # vcfdat[, P0re := fifelse(
  #   P0re > 0.95, 0.95, fifelse(P0re < 0.05, 0.05, P0re)
  # )][, P1re := fifelse(
  #   P1re > 0.95, 0.95, fifelse(P1re < 0.05, 0.05, P1re)
  # )]

  vcfdat[, logOR := log(P0 * P0re / (P1 * P1re), base = 10)]

  return(vcfdat[, hmm_fetal_hap := hmm_fetal_hap - 1])
}

#'@title getTransProbs
#'@description retrieve transition probabilities
getTransProbs <- function(map, snpset) {
  chrmap0 <- as.matrix(map[, c(2, 4)])
  chrmap1 <- cbind(snpset$pos, -1)
  chrmap1 <- chrmap1[order(chrmap1[, 1]), ]
  n0 <- nrow(chrmap0)
  n1 <- nrow(chrmap1)
  i0 <- 1
  i1 <- 1
  while (i1 <= n1 && chrmap1[i1, 1] <= chrmap0[1, 1]) {
    chrmap1[i1, 2] <- 0
    i1 <- i1 + 1
  }
  while (i1 <= n1) {
    while (i0 <= n0 && chrmap0[i0, 1] < chrmap1[i1, 1]) {
      i0 <- i0 + 1
    }
    if (i0 <= n0) {
      chrmap1[i1, 2] <- (chrmap0[i0, 2] - chrmap0[i0 - 1, 2]) *
        (chrmap1[i1, 1] - chrmap0[i0 - 1, 1]) / (chrmap0[i0, 1] - chrmap0[i0 - 1, 1]) + chrmap0[i0 - 1, 2]
      i1 <- i1 + 1
    } else {
      while (i1 <= n1) {
        chrmap1[i1, 2] <- chrmap0[n0, 2]
        i1 <- i1 + 1
      }
    }
  }
  TransProbs <- array(0, c(2, 2, n1 - 1))
  for (i in 1:(n1 - 1)) {
    p_ii1 <- (chrmap1[i + 1, 2] - chrmap1[i, 2]) / 100
    # # modifications
    # p_ii1 <- ifelse(p_ii1 < 0.0005, 0.0005, p_ii1)
    if (p_ii1 > 0.5) {
      p_ii1 <- 0.5
    }
    TransProbs[, , i] <- matrix(c(1 - p_ii1, p_ii1, p_ii1, 1 - p_ii1), 2)
  }
  return(TransProbs)
}

#'@title initHMM
#'@description initiate HMM model
initHMM <-
  function(States,
           Symbols,
           startProbs = NULL,
           transProbs = NULL,
           emissionProbs = NULL) {
    nStates <- length(States)
    nSymbols <- length(Symbols)
    S <- rep(1 / nStates, nStates)
    T <-
      0.5 * diag(nStates) + array(0.5 / (nStates), c(nStates, nStates))
    E <- array(1 / (nSymbols), c(nStates, nSymbols))
    names(S) <- States
    dimnames(E) <- list(states = States, symbols = Symbols)
    if (!is.null(startProbs)) {
      S[] = startProbs[]
    }
    if (!is.null(transProbs)) {
      if (length(dim(transProbs)) == 2) {
        T <- array(0, c(dim(transProbs), 1))
        T[, , 1] <- transProbs
        dimnames(T) <- list(from = States, to = States, NULL)
      } else if (length(dim(transProbs)) == 3) {
        T <- transProbs
        dimnames(T) <- list(from = States, to = States, NULL)
      }
    }
    if (!is.null(emissionProbs)) {
      E[, ] <- emissionProbs[, ]
    }
    return(
      list(
        States = States,
        Symbols = Symbols,
        startProbs = S,
        transProbs = T,
        emissionProbs = E
      )
    )
  }

#'@title viterbi
#'@description Viterbi algorithm for HMM decoding
viterbi <- function(hmm, observation) {
  hmm$transProbs[is.na(hmm$transProbs)] <- 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] <- 0
  nObservations <- length(observation)
  nStates <- dim(hmm$transProbs)[1]
  if (length(dim(hmm$transProbs)) == 3 &&
      dim(hmm$transProbs)[3] != nObservations - 1) {
    cat(
      "Warning:The 3rd dim of transProbs(",
      dim(hmm$transProbs)[3],
      ")
                is not as long as observations-1(",
      nObservations - 1,
      ")\n"
    )
  }

  tritransProbs <- array(rep(
    c(hmm$transProbs),
    round(1 + nStates ** 2 * nObservations / length(hmm$transProbs))
  ),
  c(nStates, nStates, nObservations - 1))
  dimnames(tritransProbs) <- dimnames(hmm$transProbs)
  nStates <- length(hmm$States)
  v <- array(NA, c(nStates, nObservations))
  dimnames(v) <-
    list(states = hmm$States, index = 1:nObservations)
  for (state in hmm$States) {
    v[state, 1] = log(hmm$startProbs[state] * hmm$emissionProbs[state, observation[1]])
  }
  if (nObservations >= 2) {
    for (k in 2:nObservations) {
      for (state in hmm$States) {
        maxi <- NULL
        for (previousState in hmm$States) {
          temp <-
            v[previousState, k - 1] + log(tritransProbs[previousState, state, k - 1])
          maxi <- max(maxi, temp)
        }
        v[state, k] <-
          log(hmm$emissionProbs[state, observation[k]]) + maxi
      }
    }
  }
  optimal_path_prob.temp <- max(v[, nObservations])
  #cat(paste0("optimal_path_prob:", optimal_path_prob.temp, "\n"))
  viterbiPath <- rep(NA, nObservations)
  for (state in hmm$States) {
    if (max(v[, nObservations]) == v[state, nObservations]) {
      viterbiPath[nObservations] <- state
      break
    }
  }
  if (nObservations >= 2) {
    for (k in (nObservations - 1):1) {
      for (state in hmm$States) {
        if (max(v[, k] + log(tritransProbs[, viterbiPath[k + 1], k])) ==
            v[state, k] + log(tritransProbs[state, viterbiPath[k + 1], k])) {
          viterbiPath[k] <- state
          break
        }
      }
    }
  }
  return(viterbiPath)
}
