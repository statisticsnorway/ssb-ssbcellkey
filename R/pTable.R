#' Probability matrix from package ptable
#'
#' The matrix of cumulative probabilities (default) can be input to \code{\link{Pconvert}}
#'
#' With doCumSum = FALSE, output is "tMatrix" created by \code{create_cnt_ptable}
#'
#' @param D \code{create_ptable} parameter
#' @param V \code{create_ptable} parameter
#' @param js \code{create_ptable} parameter
#' @param icat \code{create_ptable} parameter
#' @param table \code{create_ptable} parameter
#' @param ... \code{create_ptable} parameters
#' @param doCumSum Cumulative probabilities when TRUE
#'
#' @return Probability matrix
#' @importFrom methods slot
#' @importFrom ptable create_ptable
#' @export
#' @author Øyvind Langsrud, Daniel Lupp
#'
#' @examples
#' Pmatrix(D=3, V=3, pstay= 0.5, doCumSum = FALSE)
#' Pmatrix(D=3, V=3, pstay= 0.5)
Pmatrix <-
  function(D = 5,
           V = 3,
           js = 2,
           table = "cnts",
           icat = c(1, D),
           ...,
           doCumSum = TRUE) {
    if (!table %in% c("cnts", "nums"))
      stop("the `table` parameter must be \"cnts\" or \"nums\"")
    if (table %in% "cnts") {
      pTable <-
        ptable::create_ptable(
          D = D,
          V = V,
          js = js,
          table = table,
          ...
        )
      pMatrix <-
        slot(pTable, "tMatrix")  # slot(pTable, "pMatrix") in earlier version of ptable
      if (doCumSum)
        for (i in seq_len(nrow(pMatrix)))
          pMatrix[i, ] <- cumsum(pMatrix[i, ])
      
      if (!is.null(rownames(pMatrix)))
        if (!identical(as.integer(rownames(pMatrix)), 0:(nrow(pMatrix) - 1L)))
          warning("Unusual rownames. Pconvert may not work.")
      return(pMatrix)
    }
    pTable <-
      ptable::create_ptable(
        D = D,
        V = V,
        table = table,
        icat = icat,
        ...
      )
    return(pTable)
  }



#' Convert counts using pTable
#'
#' @param x counts or values to be perturbed
#' @param pMatrix Output from \code{\link{Pmatrix}}
#' @param rkeys  uniformly distributed keys
#' @param table character vector defining whether to perturb magnitude ("nums") or frequency("cnts") tables.'
#' @param contributors matrix  containing k top contributors, nrow(contributors) = length(x), ncol(contributors) = k.
#' @return converted counts
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' Pconvert(rep(1, 4), Pmatrix(), c(0.5, 0.8, 0.9, 0.99))
#' Pconvert(1:9, Pmatrix(), (1:9)/10)
#' Pconvert(1:9, Pmatrix(), rev((1:9)/10))
Pconvert <-
  function(x,
           pMatrix,
           rkeys = runif(NROW(x))) {
    z <- as.vector(x)
    # if (!table %in% c("cnts", "nums"))
    #   stop("the `table` parameter must be \"cnts\" or \"nums\"")
    # if (table == "cnts") {
    z[z >= nrow(pMatrix)] <- nrow(pMatrix) - 1L
    pM1 <- pMatrix[z + 1, , drop = FALSE]
    cc <- col(pM1) * as.integer(!(pM1 < rkeys))
    cc[cc == 0] <-
      .Machine$integer.max # Instead of Inf to keep integer class
    return(x - z + as.integer(colnames(pMatrix))[apply(cc, 1, min)])
    # } else {
    #   if (is.null(contributors))
    #     stop("contributors matrix cannot be NULL for table = \"nums\"")
    #   a0 <- as.numeric(rownames(pMatrix)[1])
    #   a1 <- as.numeric(rownames(pMatrix)[nrow(pMatrix)])
    #   return(c(a0, a1))
    # }
  }
#
# contributors <- matrix(6:1, ncol = 3)
# num <- rowSums(contributors) + c(4,2)
# pm <- Pmatrix(table = "nums", step = 2)
# m <- c(0.8,0.7,0.6)
# Pconvert_num <- function(num, pMatrix, m, contributors, rkeys = matrix(runif(ncol(contributors) * nrow(contributors)), ncol = ncol(contributors))) {
#   # return(contributors / contributors %*% diag(m))
#   lambda <- contributors %*% diag(1/m)
#   D <- max(as.integer(rownames(pMatrix)))
#   lambda[] <- vapply(lambda, function(x) min(x, D), numeric(1))
#   # here we set a0 to 1, a1 to D
#   lambda <- (lambda - 1) / (D-1)
#
#   pM1 <- pMatrix[c(1, D) + 1, , drop = FALSE]
#   cc <- col(pM1) * as.integer(!(pM1 < rkeys))
#   cc[cc == 0] <-
#     .Machine$integer.max
#   as.numeric(colnames(pMatrix))[apply(cc, 1, min)]
# }
# l <- Pconvert_num(num = num, pMatrix = pm, m = m, contributors =  contributors)
# l
#' pMatrix
