
#' CellKeyFits
#'
#' @param data data frame (inner cells)
#' @param freqVar Variable holding counts 
#' @param rKeyVar Variable holding keys 
#' @param formula Model formula  
#' @param hierarchies List of hierarchies
#' @param dimVar Dimensional variables
#' @param preAggregate Aggregation
#' @param extend0  Extend0  
#' @param iter Mipf parameter
#' @param eps  Mipf parameter
#' @param tol  Mipf parameter
#' @param reduceBy0 Mipf parameter
#' @param reduceByColSums Mipf parameter
#' @param reduceByLeverage Mipf parameter
#' @param ... dots
#' 
#' @return
#' @importFrom SSBtools Mipf
#' @importFrom Matrix crossprod
#' @export
#'
#' @examples
#' 1 + 1 
CellKeyFits <- function(data, freqVar = NULL, rKeyVar = NULL, hierarchies = NULL, formula = NULL, dimVar = NULL, 
                        preAggregate = is.null(freqVar), extend0 = TRUE, iter = 1000, eps = 0.01,
                        tol = 1e-13, reduceBy0 = TRUE, reduceByColSums = TRUE, reduceByLeverage = FALSE, ...) {
  
  
  freqVar <- names(data[1, freqVar, drop = FALSE])
  
  
  if (preAggregate | extend0) {
    
    # 1
    data <- CellKey(data = data, freqVar = freqVar, rKeyVar = rKeyVar, 
                    hierarchies = hierarchies, formula = formula, dimVar = dimVar, 
                    preAggregate = preAggregate, innerReturn = 1, ...)
    
    
    ma <- match(c("freqVar", "f_Re_qVa_r"), names(data))
    ma <- ma[!is.na(ma)]
    freqVar <- c(names(data)[ma], "freq")[1]
    
    nrowOrig <- nrow(data)
    
    # 2
    if (extend0) {
      data <- Extend0(data, freqVar)
    }
  } else {
    nrowOrig <- nrow(data)
  }
  
  
  # 3
  mm <- ModelMatrix(data = data, hierarchies = hierarchies, formula = formula, crossTable = TRUE, dimVar = dimVar, ...)
  
  
  # 4
  if (nrowOrig < nrow(data)) {
    a <- CellKey(data = data[seq_len(nrowOrig), ], freqVar = freqVar, rKeyVar = rKeyVar, preAggregate = FALSE, 
                 x = mm$modelMatrix[seq_len(nrowOrig), , drop = FALSE], crossTable = mm$crossTable, ...)
  } else {
    a <- CellKey(data = data, freqVar = freqVar, rKeyVar = rKeyVar, preAggregate = FALSE, 
                 x = mm$modelMatrix, crossTable = mm$crossTable, ...)
  }
  
  
  # 5
  lsFit <- LSfitNonNeg(mm$modelMatrix, Matrix(a$perturbed, ncol = 1))
  
  
  # 6
  ipFit <- Mipf(mm$modelMatrix, z = lsFit, iter = iter, eps = eps, tol = tol, 
                reduceBy0 = reduceBy0, reduceByColSums = reduceByColSums, reduceByLeverage = reduceByLeverage)
  
  
  list(inner = cbind(data, ipFit = as.vector(ipFit)), 
       publish = cbind(a, lsFit = as.vector(lsFit), ipFit = as.vector(crossprod(mm$modelMatrix, ipFit))))
  
}                           
 




