#' generate a random population
#'
#' @param n census size (default: 100)
#' @param ploidy ploidy of all samples (default: 1)
#' @param freq an optional vector of named allele frequencies.
#' @param nall an optional vector of allele numbers per locus
#' @param clone_gen number of generations to undergo clonal reproduction
#' @param mu mutation rate (number of mutations per generation over the whole
#'   population)
#' @param verbose when `TRUE` (default), messages will be displayed as clonal
#'   simulation progresses.
#' @param genclone when `TRUE` (default), data will be converted to genclone.
#'
#' @return a `genind` or `genclone` object.
#' @export
#'
#' @details By default, this will simulate 11 loci from 4 to 12 alleles/locus
#'   using allele frequencies drawn from a uniform distribution. Populations are
#'   initially created from a multinomial distribution at each locus separately.
#'
#'   Clonal reproduction consists of first mutating a single individual at a
#'   single locus by shuffling the alleles at that locus. After mutation, the
#'   population is randomly sampled with replacement.
#'
#' @examples
#' pop_generator() # panmictic population
#' pop_genearator(clone_gen = 5, mu = 0.5)
#' 
#' @importFrom poppr as.genclone
#' @importFrom adegenet genind
#' @importFrom purrr map_chr map
#' @importFrom dplyr progress_estimated
pop_generator <- function(n = 1000, 
                          samp = 100,
                          ploidy = 1, 
                          freq = NULL,
                          nall = sample(4:12, 11, replace = TRUE),
                          mate_gen = 100,
                          clone_gen = 0,
                          mu = 0.05,
                          verbose = TRUE,
                          genclone = TRUE){
  if (!is.null(freq)) {
    loc_all <- strsplit(names(freq), "\\.")
    locnames <- purrr::map_chr(loc_all, 1)
    alleles  <- purrr::map_chr(loc_all, 2) %>% as.integer()
    names(freq) <- alleles
    freqlist <- split(freq, locnames)
  } else {
    freqlist <- purrr::map(nall, ~{setNames(runif(.), seq(.))})
    freqlist <- purrr::map(freqlist, ~{./sum(.)})
    names(freqlist) <- paste("locus", seq(nall))
    locnames <- rep(names(freqlist), lengths(freqlist))
  }
  reslist <- purrr::map(freqlist, ~{ t(rmultinom(n, ploidy, .)) })
  names(reslist) <- names(freqlist)
  for (i in names(reslist)) {
    colnames(reslist[[i]]) <- paste(i, colnames(reslist[[i]]), sep = ".")
  }
  restab <- do.call(cbind, reslist)
  mutation_events <- 0
  locfac <- factor(locnames, unique(locnames))
  if (verbose) {
    message("Beginning mating\n")
    p <- dplyr::progress_estimated(mate_gen)
  } 
  for (i in seq_len(mate_gen)) {
    if (verbose) p$tick()$print()
    if (mu > runif(1)) {
      restab <- mutator(restab, locnames, ploidy)
      mutation_events <- mutation_events + 1
    }
    restab <- t(vapply(seq_len(n), 
                       FUN = function(i){ mater(restab, locfac, ploidy) }, 
                       FUN.VALUE = restab[1, ]))
  }
  if (verbose) {
    message("\nI recorded ", mutation_events, " mutation events\n")
    p$stop()
  }
  
  if (clone_gen > 0) {
    if (verbose) {
      message("starting clonal reproduction...\n")
      p <- dplyr::progress_estimated(clone_gen)
    }
    for (i in seq_len(clone_gen)) {
      if (verbose) p$tick()$print()
      if (mu > runif(1)) {
        restab <- mutator(restab, locnames, ploidy)
        mutation_events <- mutation_events + 1
      }
      restab <- restab[sample(n, replace = TRUE), ]
    }
    if (verbose) {
      message("\nI recorded ", mutation_events, " mutation events\n")
      p$stop()
    }
  }
  res <- genind(restab[sample(n, samp), ], ploidy = ploidy, type = "codom")
  if (genclone) {
    res <- as.genclone(res)
  }
  res
}

#' Create a single recombinant offspring from two randomly sampled parents
#'
#' @param mat a matrix of individuals x alleles
#' @param loci a factor deliniating the loci in each column
#' @param ploidy the ploidy of the sample
#'
#' @return a vector of alleles
#' @export
#'
#' @examples
#' mat <- t(rbind(rmultinom(100, 2, runif(10)),
#'              rmultinom(100, 2, runif(10)))) # two loci with 20 alleles
#' loc <- factor(rep(LETTERS[1:2], each = 10))
#' t(vapply(seq_len(nrow(mat)), function(i) mater(mat, loc, 2), mat[1, ]))
mater <- function(mat, loci, ploidy = 1){
  x        <- sample(nrow(mat), 2)
  parentab <- colSums(mat[x, ])
  res      <- rep(0L, ncol(mat))
  for (i in seq_len(ploidy)) {
    pindex   <- which(parentab > 0)
    pl       <- split(parentab[pindex], loci[pindex])
    new_nall <- lengths(pl)
    samps    <- vapply(new_nall, sample, integer(1), 1) + cumsum(new_nall) - new_nall
    res[pindex[samps]] <- 1L
    parentab[pindex[samps]] <- parentab[pindex[samps]] - 1L 
  }
  res
}

mutator <- function(mat, locnames, ploidy){
  ind <- sample(nrow(mat), 1)
  loc <- grepl(sample(locnames, 1), colnames(mat))
  prob <- colMeans(mat[, loc], na.rm = TRUE)
  mat[ind, loc] <- rmultinom(1, ploidy, 1 - prob)
  mat
}
