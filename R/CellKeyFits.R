
#' CellKeyFits
#'
#' @param data data frame (inner cells)
#' @param freqVar Variable holding counts 
#' @param rKeyVar Variable holding keys 
#' @param formula Model formula  
#' @param hierarchies List of hierarchies
#' @param dimVar Dimensional variables
#' @param preAggregate Aggregation
#' @param xReturn	Dummy matrix in output when `TRUE`. To return crossTable as well, use `xReturn = 2`.
#' @param extend0  Extend0
#' @param limit LSfitNonNeg parameter
#' @param viaQR LSfitNonNeg parameter 
#' @param iter Mipf parameter
#' @param eps  Mipf parameter
#' @param tol  Mipf parameter
#' @param reduceBy0 Mipf parameter
#' @param reduceByColSums Mipf parameter
#' @param reduceByLeverage Mipf parameter
#' @param ... Further parameters to CellKey
#' 
#' @return list
#' @importFrom SSBtools Mipf Extend0 LSfitNonNeg DummyDuplicated
#' @importFrom Matrix crossprod
#' @export
#'
#' @examples
#' my_km2 <- SSBtools::SSBtoolsData("my_km2")
#' CellKeyFits(my_km2, "freq", formula = ~(Sex + Age) * Municipality * Square1000m + Square250m)  
CellKeyFits <- function(data, freqVar = NULL, rKeyVar = NULL, hierarchies = NULL, formula = NULL, dimVar = NULL, 
                        preAggregate = is.null(freqVar), xReturn = FALSE, 
                        extend0 = TRUE, 
                        limit = 1e-10, viaQR = FALSE,
                        iter = 1000, eps = 0.01,
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
      extraVar <- c(rKeyVar, "rKeyVar")[1]
      if (!(extraVar %in% names(data))) {
        extraVar <- FALSE
      } 
      data <- Extend0(data, freqVar,  extraVar = extraVar)
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
  
  dd <- DummyDuplicated(mm$modelMatrix)
  
  # 5
  lsFit <- LSfitNonNeg(mm$modelMatrix[, !dd, drop=FALSE], Matrix(a$perturbed[!dd], ncol = 1), limit = limit, viaQR = viaQR)
  
  
  # 6
  ipFit <- Mipf(mm$modelMatrix, z = lsFit, iter = iter, eps = eps, tol = tol, 
                reduceBy0 = reduceBy0, reduceByColSums = reduceByColSums, reduceByLeverage = reduceByLeverage)
  
  
  data <- list(inner = cbind(data, ipFit = as.vector(ipFit)), 
               publish = cbind(a, lsFit = as.vector(lsFit), ipFit = as.vector(crossprod(mm$modelMatrix, ipFit))))
  
  if (xReturn) {
    names(mm)[1] <- "x"
    if (xReturn != 2) {
      mm <- mm[1]
    }
    return(c(data, mm))
  }
  
  data
}                           
 




