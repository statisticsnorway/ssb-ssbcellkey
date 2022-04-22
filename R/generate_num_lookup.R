#' Function for generating a probability matrix for cell key perturbation of
#' numeric variables.
#'
#'
#' @param D numeric scalar, describing maximum noise
#' @param step numeric scalar, determinng step width for noise
#' @param dcomponent data driven component used to calculate noise distribution
#' @param noisefunction function to calculate noise distribution among elements
#'  of noisefactor, based on dcomponent.
#' @param minpositive numeric scalar describing lowest positive value after
#' perturbation. Only tested with frequencies and maxentropy.
#' @param percnoise boolean vector of length 1. If TRUE, noise and data
#' components are divided by 100 in the output in order to represent percentages
#'  (technical fix to ensure no bugs because of machine precision).
#' @param ...
#'
#' @return Output is a mxn matrix, where rows correspond to m possible noise
#' factors, and columns correspond to n possible data-driven components.
#' @export
#'
#' @author Daniel P. Lupp, Hege Marie Bøvelstad
#'
#' @examples
#' a1 <- generate_prob_matrix(D = 3,
#'                           step = 1,
#'                           dcomponent = 0:5,
#'                           noisefunction = normal_noise)
#' a2 <- generate_prob_matrix(D = 3,
#'                           step = 1,
#'                           dcomponent = 1:5,
#'                           noisefunction = freq_maxentropy)
#' a3 <- generate_prob_matrix(D = 3,
#'                           step = 1,
#'                           dcomponent = 1:5,
#'                           minpositive = 2,
#'                           noisefunction = freq_maxentropy)
#'
#' dcomponent <- seq(from = 50, to = 100, by = 10)
#' a4 <-  generate_prob_matrix(
#'     D = 15,
#'     step = 1,
#'     dcomponent,
#'     split_triangular,
#'     width = 5,
#'     ddc2noise = function(x)
#'       x / 10,
#'     percnoise = TRUE
#'   )
#'
#' Visualization of noise distributions, feel free to vary column index in each
#' example
#' plot(x = as.numeric(rownames(a1)), y = a1[,1])
#' plot(x = as.numeric(rownames(a2)), y = a2[,4], type = "l")
#' plot(x = as.numeric(rownames(a3)), y = a3[,1])
#' plot(x = as.numeric(rownames(a4)), y = a4[,4])
generate_prob_matrix <- function(D,
                                 step,
                                 dcomponent,
                                 noisefunction,
                                 minpositive = 1,
                                 percnoise = FALSE,
                                 ...) {
  if (!is.function(noisefunction))
    stop("The parameter noisefunction must be a function.")
  if (!("..." %in% names(formals(noisefunction))))
    stop("You must add \"...\" to the noise function parameter list.")
  noisefactor <- seq(from = -D, to = D, by = step)
  # prob_table is output matrix
  prob_matrix <-
    matrix(NA,
           nrow = length(noisefactor),
           ncol = length(dcomponent))
  rownames(prob_matrix) <- noisefactor
  colnames(prob_matrix) <- dcomponent
  for (i in dcomponent) {
    prob_matrix[, match(i, dcomponent)] <-
      noisefunction(ddc = i, noisefactor = noisefactor, minpositive, ...)
  }
  if (percnoise) {
    rownames(prob_matrix) <- noisefactor / 100
    colnames(prob_matrix) <- dcomponent / 100
  }
  prob_matrix
}


#' Generate a lookup table given noise probabilities.
#'
#' @param probmatrix matrix containing noise probabilities
#' @param ncellkeys number of cellkeys.
#' @param datacomponents data-driven components used for lookup
#' @param seed set seed for reproducibility of random results.
#'
#' @return returns a mxn matrix, for m possible cell keys and n possible data
#' driven components used for lookup in the cell key method.
#' @export
#'
#' @author Daniel P. Lupp, Hege Marie Bøvelstad
#'
#' @examples
#' dcomponent <- seq(from = 50, to = 100, by = 10)
#' a4 <-  generate_prob_matrix(
#'     D = 15,
#'     step = 1,
#'     dcomponent,
#'     split_triangular,
#'     width = 3,
#'     ddc2noise = function(x)
#'       x / 10,
#'     percnoise = TRUE
#'   )
#' a <- generate_lookup_table(a4, 16, dcomponent/100)
generate_lookup_table <- function(probmatrix,
                                  ncellkeys,
                                  datacomponents,
                                  seed = 123) {
  set.seed(seed)
  x <- as.numeric(rownames(probmatrix))
  nc <- ncol(probmatrix)
  m <- 0
  lookup <- NULL
  while ((m <- m + 1) <= nc) {
    prob <- probmatrix[, m]
    z <- sample(x, ncellkeys, TRUE, prob)
    lookup <- cbind(lookup, z)
  }
  m <- nc
  prob <-  probmatrix[, nc]
  ndatacomponents <- length(datacomponents)
  while ((m <- m + 1) <= ndatacomponents) {
    z <- sample(x, ncellkeys, TRUE, prob)
    lookup <- cbind(lookup, z)
  }
  colnames(lookup) <- datacomponents
  rownames(lookup) <- 0:(ncellkeys - 1)
  lookup
}

retrieve_noise <- function(prob_table, keys, ddc) {
  key_lookup <- match(keys, rownames(prob_table))
  ddc_lookup <- match(ddc, colnames(prob_table))
  mapply(function(x, y) prob_table[x, y],
         key_lookup, ddc_lookup)
}

perturb <- function(keys, values, ptable, ddc.function = function(x) x) {
  ddc <- ddc.function(values)
  values + retrieve_noise(ptable, keys, ddc)
}

truncate_int_by_bit <- function(ints, nobits = 8) {
  if (!length(ints))
    return(NULL)
  b <- sapply(ints, function(x)
    as.integer(intToBits(x)))
  nobits <- min(32, nobits)
  b <- matrix(b[1:nobits,], ncol = length(ints))
  apply(b, 2, function(y)
    ifelse(all(y == 0), 0, Reduce(sum, sapply(which(y > 0) - 1, function(x)
      2 ^ x))))
}

aggregate_bitkeys <- function(rkeys,
                             pop.key = NULL,
                             trunc.bit = 8) {
  if (!is.null(pop.key))
    return(bitwXor(pop.key,
                   as.integer(Reduce(
                     bitwXor,
                     truncate_int_by_bit(rkeys, nobits = trunc.bit)
                   ))))
  Reduce(bitwXor, truncate_int_by_bit(rkeys, nobits = trunc.bit))
}
