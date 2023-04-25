#' Probability matrix from package ptable  
#' 
#' The matrix of cumulative probabilities (default) can be input to \code{\link{Pconvert}}
#' 
#' With doCumSum = FALSE, output is "pMatrix" created by \code{create_cnt_ptable}
#'
#' @param D \code{create_cnt_ptable} parameter
#' @param V \code{create_cnt_ptable} parameter
#' @param js \code{create_cnt_ptable} parameter
#' @param ... \code{create_cnt_ptable} parameters
#' @param doCumSum Cumulative probabilities when TRUE
#'
#' @return Probability matrix
#' @importFrom methods slot
#' @importFrom ptable create_cnt_ptable
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' Pmatrix(D=3, V=3, pstay= 0.5, doCumSum = FALSE)
#' Pmatrix(D=3, V=3, pstay= 0.5)
Pmatrix <- function(D = 5, V = 3, js = 2,..., doCumSum = TRUE) {
  if (requireNamespace("ptable", quietly = TRUE)) {
    pTable <- ptable::create_cnt_ptable(D = D, V = V, js = js, ...)
    pMatrix <- slot(pTable, "pMatrix")
    if (doCumSum) 
      for (i in seq_len(nrow(pMatrix))) pMatrix[i, ] <- cumsum(pMatrix[i, ])
  } else {
    pMatrix <- matrix(c(1, 0.37, 0.2, 1, 0.73, 0.4, 1, 0.9, 0.6, 1, 1, 0.8, 1, 1, 1), 3, 5, dimnames = list(0:2, 0:4))
    warning("Package ptable not available. \"round(Pmatrix(D=2,js=0),2)\" returned. Input parameters ignored.")
  }
  if (!is.null(rownames(pMatrix))) 
    if (!identical(as.integer(rownames(pMatrix)), 0:(nrow(pMatrix) - 1L))) 
      warning("Unusual rownames. Pconvert may not work.")
  pMatrix
}



#' Convert counts using pTable
#'
#' @param x counts
#' @param pMatrix Output from \code{\link{Pmatrix}}
#' @param rkeys  uniformly distributed keys
#'
#' @return converted counts
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' Pconvert(rep(1, 4), Pmatrix(), c(0.5, 0.8, 0.9, 0.99))
#' Pconvert(1:9, Pmatrix(), (1:9)/10)
#' Pconvert(1:9, Pmatrix(), rev((1:9)/10))
Pconvert <- function(x, pMatrix, rkeys = runif(NROW(x))) {
  z <- as.vector(x)
  z[z >= nrow(pMatrix)] <- nrow(pMatrix) - 1L
  pM1 <- pMatrix[z + 1, , drop = FALSE]
  cc <- col(pM1) * as.integer(!(pM1 < rkeys))
  cc[cc == 0] <- .Machine$integer.max # Instead of Inf to keep integer class 
  x - z + as.integer(colnames(pMatrix))[apply(cc, 1, min)]
}


