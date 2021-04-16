## Adapted from https://github.com/dieterich-lab/DesignMetabolicRNAlabeling


##' Kinetic definitions for the biochemical separation 
##' Pulse model
makeFormStdPulse <- function(fractions, contaminated=FALSE) {
  formulas <- MeanFormulas(
    total      = exp(mu),
    labelled   = exp(mu) * (1 - exp(-d * time)),
    unlabelled = exp(mu - d * time)
  )
  formulaIndexes  <- list(
    total         = "total",
    flow_through  = c("unlabelled", "labelled"),
    pull_down     = c("labelled", "unlabelled")
  )[fractions]
  if (!contaminated) {
    formulaIndexes <- lapply(formulaIndexes, `[[`, 1)
  }
  usedFormulas <- unique(unlist(formulaIndexes))
  list(formulas    = formulas[usedFormulas],
       formulaIndexes = formulaIndexes)
}


##' Kinetic definitions for nucleotide conversion  
##' Pulse model
##' mu1 = background unlabelled
##' mu2 = background labelled
##' mu3 = difference between maximum and background labelled
##' The total fraction is mu = mu1 + mu2 + mu3
##' with mu1 > mu2 + mu3 and mu3 > 0 (lower bound)
##' unlabelled = mu1 + mu3 * exp(-d * time)
##' labelled = mu2 + mu3 * (1 - exp(-d * time))
##' Reparametrizing yields:
makeFormRtcPulse <- function() {
  list(
    formulas = MeanFormulas(
      unlabelled = exp(mu1) + exp(mu2) + exp(mu3) * (1 + exp(-d * time)),
      labelled = exp(mu2) + exp(mu3) * (1 - exp(-d * time))),
    formulaIndexes = list(
      unlabelled = "unlabelled",
      labelled = "labelled"))
}

