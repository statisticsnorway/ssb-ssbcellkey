
#' Converting tabular data by a perturbation table and uniformly distributed keys
#' 
#' The perturbation table is created by package ptable. 
#' 
#' With `freqVar=NULL` and `preAggregate = FALSE`, the calculated keys are returned.
#' 
#' @param data Input data as a data frame (inner cells)
#' @param freqVar Variable holding counts (name or number)
#' @param rKeyVar Variable holding uniformly distributed keys
#' @param formula Model formula defining publishable cells. Will be used to calculate \code{x} (via \code{\link{ModelMatrix}}). 
#' @param hierarchies List of hierarchies, which can be converted by \code{\link{AutoHierarchies}}. 
#'        Thus, a single string as hierarchy input is assumed to be a total code. 
#'        Exceptions are \code{"rowFactor"} or \code{""}, which correspond to only using the categories in the data.
#' @param dimVar The main dimensional variables and additional aggregating variables. This parameter can be  useful when hierarchies and formula are unspecified. 
#' @param preAggregate When `TRUE`, the data will be aggregated within the function to an appropriate level. 
#' @param total	String used to name totals
#' @param x Dummy matrix defining cells to be published (possible as input instead of generated)
#' @param crossTable	Data frame to accompany `x` when `x` is input.  
#' @param D \code{\link{pt_create_pParams}} parameter
#' @param V \code{\link{pt_create_pParams}} parameter
#' @param js \code{\link{pt_create_pParams}} parameter
#' @param pstay \code{\link{pt_create_pParams}} parameter
#'
#' @return Data frame with keys or aggregated counts (original and perturbated) 
#' 
#' 
#' @importFrom SSBtools FindHierarchies ModelMatrix
#' @importFrom stats runif aggregate as.formula delete.response terms
#' @importFrom utils flush.console
#' @importFrom Matrix Matrix   
#' @export
#'
#' @examples
#' z <- data.frame(geo  = c("Iceland", "Portugal", "Spain"), 
#'                 eu = c("nonEU", "EU", "EU"),
#'                 year = rep(c("2018","2019"), each = 3),
#'                 freq = c(2,3,7,1,5,6), stringsAsFactors = FALSE)
#' 
#' CellKeyViaDummy(z, "freq", formula = ~eu * year + geo)
#' CellKeyViaDummy(z, freqVar = NULL, formula = ~eu * year + geo)
#' CellKeyViaDummy(z, freqVar = NULL, formula = ~eu * year + geo, preAggregate = FALSE)
#' z$keys <- sin(1:6)%%1
#' CellKeyViaDummy(z, "freq", "keys", formula = ~eu * year + geo)
#' CellKeyViaDummy(z, "freq", "keys", dimVar = c("geo", "eu", "year"))
#' CellKeyViaDummy(z, freqVar = NULL, "keys", dimVar = c("geo", "eu", "year"))
#' CellKeyViaDummy(z, "freq", "keys", hierarchies = 
#'      list(geo = c("EU", "@Portugal", "@Spain", "Iceland"), year = c("2018", "2019")))
#'      
CellKeyViaDummy <- function(data, freqVar=NULL, rKeyVar=NULL, 
                            hierarchies = NULL, formula = NULL, dimVar = NULL, 
                            preAggregate = is.null(freqVar),
                            total = "Total", 
                            x = NULL, crossTable = NULL,
                            D=5, V=3, js=2, pstay = NULL){

  force(preAggregate)
    
  # Ensure character (integer possible input) 
  freqVar <- names(data[1, freqVar, drop = FALSE])
  
  if (is.null(formula) & is.null(hierarchies) & is.null(x) & is.null(dimVar)){
    dimVar <- names(data[1, !(names(data) %in% c(freqVar, rKeyVar)), drop = FALSE])
  } else {
    dimVar <- names(data[1, dimVar, drop = FALSE])
  }
  
  
  if (preAggregate) {
    cat("[preAggregate ", dim(data)[1], "*", dim(data)[2], "->", sep = "")
    flush.console()
    if (!is.null(hierarchies)) {
      dVar <- names(hierarchies)
    } else {
      if (!is.null(formula)) {
        dVar <- row.names(attr(delete.response(terms(as.formula(formula))), "factors"))
      } else {
        dVar <- dimVar
      }
    }
    if(length(c(freqVar, rKeyVar))){
      if(!length(freqVar)){  # simpler code by adding 1s, but then the entire data.frame is copied into memory
        freqVar_ <- "f_Re_qVa_r"
        data_freqVar <- aggregate(list(f_Re_qVa_r = data[[dVar[1]]]), data[, dVar, drop = FALSE], length)[[freqVar_]]
      } else {
        data_freqVar <- NULL
      }
      data <- aggregate(data[, c(freqVar, rKeyVar), drop = FALSE], data[, dVar, drop = FALSE], FUNaggregate)
      if(!is.null(data_freqVar)){
        freqVar <- freqVar_
        data[[freqVar]] <- data_freqVar
      }
    } else {
      data <- aggregate(list(f_Re_qVa_r = data[[dVar[1]]]), data[, dVar, drop = FALSE], length)
      freqVar <- "f_Re_qVa_r" 
    }
    cat(dim(data)[1], "*", dim(data)[2], "]", sep = "")
    flush.console()
  }
  if (is.null(rKeyVar)) {
    rkeys <- runif(NROW(data))
    rKeyVar <- "rKeyVar"
  } else {
    rkeys <- data[, rKeyVar]
  }
  
  
  if (is.null(x)) {
    if (!is.null(crossTable)) {
      warning("crossTable in input ignored when input x is NULL")
    }
    crossTable <- TRUE
  }
  
  cat("[ModelMatrix")
  flush.console()
  
  
  mm <- ModelMatrix(data = data, hierarchies = hierarchies, formula = formula, crossTable = crossTable, modelMatrix = x, total = total, dimVar = dimVar)
  
  
  cat("]")
  flush.console()
  
  cat("[aggregate rkeys")
  flush.console()
  # rKey <- Matrix::crossprod(mm$modelMatrix, rkeys)[, 1, drop = TRUE] %%1 
  # More general calculation with DummyApply and FUNaggregate
  rKey <- DummyApply(mm$modelMatrix, rkeys, FUNaggregate)
  cat("]")
  flush.console()
  
  
  cat("[Aggregates.")
  flush.console()
  
  if (length(freqVar)) {
    yOrig <- Matrix::crossprod(mm$modelMatrix, data[, freqVar])[, 1, drop = TRUE]
  } else {
    data <- cbind(as.data.frame(mm$crossTable, stringsAsFactors = FALSE), r_Ke_yVa_r = rKey)
    names(data)[NCOL(data)] <- rKeyVar
    cat("]\n")
    flush.console()
    rownames(data) <- NULL
    return(data)
  }
  
  cat(".ptable.")
  flush.console()
  
  if (is.character(rKey)) {
    bitkey <- matrix(rKey, ncol = 1, dimnames = list(NULL, rKeyVar))
    rKey <- (as.numeric(rKey) + 0.5)/256
    rKeyVar <- paste0("unif_", rKeyVar)
  } else {
    bitkey <- matrix("0", nrow = length(rKey), ncol = 0)
  }
  
  pMatrix <- Pmatrix(D = D, V = V, js = js, pstay = pstay)
  
  cat(".")
  flush.console()
  
  yPert <- Pconvert(yOrig, pMatrix, rKey)
  
  cat(".")
  flush.console()
  
  data <- cbind(as.data.frame(mm$crossTable, stringsAsFactors = FALSE), yOrig = yOrig, yPert = yPert, bitkey, r_Ke_yVa_r = rKey)
  cat("]\n")
  flush.console()
  names(data)[NCOL(data)] <- rKeyVar
  rownames(data) <- NULL
  data
}




FUNaggregate <- function(x){ 
  if(is.character(x)){ 
    stop("Character input is not implemented yet.")
  }
  if(!length(x))
    return(sum(x))   ### OBS:  No rKey -> 0
  if(x[1]>0 & x[1]<1)
    return(sum(x) %%1)
  sum(x)
}




















