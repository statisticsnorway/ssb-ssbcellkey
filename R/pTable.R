#' Probability matrix from package ptable  
#' 
#' The matrix of cumulative probabilities (default) can be input to \code{\link{Pconvert}}
#' 
#' With doCumSum = FALSE, output is "tMatrix" created by \code{create_cnt_ptable}
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
    pTable <- ptable::create_cnt_ptable(D = D, V = V, js = js, ...)
    pMatrix <- slot(pTable, "tMatrix")  # slot(pTable, "pMatrix") in earlier version of ptable 
    if (doCumSum) 
      for (i in seq_len(nrow(pMatrix))) pMatrix[i, ] <- cumsum(pMatrix[i, ])
  
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


