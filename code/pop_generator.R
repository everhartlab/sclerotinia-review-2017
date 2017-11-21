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
pop_generator <- function(n = 100, 
                          ploidy = 1, 
                          freq = NULL,
                          nall = sample(4:12, 11, replace = TRUE),
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
  if (clone_gen > 0) {
    if (verbose) p <- dplyr::progress_estimated(clone_gen)
    for (i in seq_len(clone_gen)) {
      if (verbose) p$tick()$print()
      if (mu > runif(1)) {
        ind <- sample(n, 1)
        loc <- grepl(sample(locnames, 1), colnames(restab))
        prob <- colMeans(restab[, loc], na.rm = TRUE)
        restab[ind, loc] <- rmultinom(1, ploidy, prob)
        mutation_events <- mutation_events + 1
      }
      restab <- restab[sample(n, replace = TRUE), ]
    }
    if (verbose) {
      cat("I recorded", mutation_events, "mutation events\n")
      p$stop()
    }
  }
  res <- genind(restab, ploidy = ploidy, type = "codom")
  if (genclone) {
    res <- as.genclone(res)
  }
  res
}

message("pop_generator is loaded.\n\nUSAGE:\n")
x <- capture.output(args(pop_generator))
x <- gsub("NULL", "", x)
x <- gsub("function ", "pop_generator", x)
x <- gsub(" ,", " NULL,", x)
message(paste(trimws(x), collapse = " "))
