#' Uniform noise function for perturbation using the cell-key method 
#'
#' @param noisefactor 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
uniform_noise <- function(noisefactor,
                          ...) {
  output <- rep(1 / length(noisefactor), length(noisefactor))
  names(output) <- noisefactor
  output
}

#' Normal noise function for perturbation using the cell-key method
#'
#' @param noisefactor 
#' @param ddc 
#' @param step 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
normal_noise <- function(noisefactor, 
                          ddc,
                          step = step,
                          ...){
  output <-
    truncnorm::dtruncnorm(
      noisefactor,
      a = min(noisefactor),
      b = max(noisefactor),
      mean = 0,
      sd = 1
    ) * step
  names(output) <- noisefactor
  output
}

laplace_noise <- function(noisefactor,
                          ddc,
                          split = 0,
                          scale = 2,
                          ...) {
     if (is.function(split))
      split <- sapply(ddc, split)
  output <- rep(0, length(noisefactor))
  names(output) <- noisefactor
  m <- max(abs(max(noisefactor)), abs(min(noisefactor)))
  m <- m-split
  nf <- noisefactor[abs(noisefactor) <= m]
  l <- LaplacesDemon::dtrunc(nf, 
                        "laplace",
                        a = min(nf), 
                        b = max(nf),
                        location = 0,
                        scale = scale)
  output[noisefactor <= -split] <- l[1:(floor(length(l)/2) + 1)]
  output[noisefactor >= split] <- l[(floor(length(l)/2) + 1):length(l)]
  output
}

#' Symmetric two peak triangular noise.
#'
#' @param noisefactor vector of possible noise factors.
#' @param ddc data driven component used for lookup.
#' @param width numeric scalar, width of noise triangle.
#' @param ddc2noise function to translate `ddc` to the distance between
#' triangles and 0. If set to NULL (default), `ddc` is equal to distance between
#' triangles and 0.
#' @param ...
#'
#' @return numeric vector describing noise distribution.
#' @export
#' 
#'
#' @examples
#' nf <- seq(from = -15, to = 15, by = 1)
#' p1 <- split_triangular(noisefactor = nf,
#'                        ddc = 5,
#'                        width = 4)
#' p2 <- split_triangular(noisefactor = nf,
#'                        ddc = 5,
#'                        ddc2noise = function(x) 2*x,
#'                        width = 4)
#' plot(x = nf, y = p1)
#' plot(x = nf, y = p2)
split_triangular <- function(noisefactor,
                             ddc,
                             width = ddc/3,
                             ddc2noise = NULL,
                             ...) {
  if (!is.null(ddc2noise) & is.function(ddc2noise))
    ddc <- sapply(ddc, ddc2noise)
  if (abs(ddc) > max(abs(noisefactor)))
    stop(
      "The data component is not aligned with possible noise factors. Please reconsider your parameters."
    )
  if (ddc == max(noisefactor))
    stop(
      "The data component is the same as max noise factor, no width is possible. Please reconsider your parameters."
    )
  if (ddc > 0 & ddc + width > max(noisefactor)) {
    warning("Too wide noise, adjusted to fit to noise factors.")
    width <- max(noisefactor) - ddc
  }
  noisefactor <- sort(noisefactor)
  output <- rep(0, length(noisefactor))
  negs <-
    which(noisefactor >= (-ddc - width) & noisefactor <= -ddc)
  pos <- which(noisefactor <= (ddc + width) & noisefactor >= ddc)
  length_negs <- length(negs)
  length_pos <- length(pos)
  output[negs] <-
    sapply(seq_len(length_negs), function(x)
      x / ((length_negs) * (length_negs + 1)))
  output[pos] <-
    rev(sapply(seq_len(length_pos), function(x)
      x / ((length_negs) * (length_negs + 1))))
  output
}

#' Maximum entropy noise function for frequency table perturbation using the 
#' cell-key method
#'
#' @param noisefactor 
#' @param ddc 
#' @param minpositive 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
freq_maxentropy <- function(noisefactor,
                            ddc,
                            minpositive,
                            ...) {
  prob <- calcprob(b = minpositive,
                   d = max(noisefactor),
                   n = ddc)
  output <- rep(0, length(noisefactor))
  names(output) <- noisefactor
  output[match(prob$X, names(output))] <- prob$p
  output
}

#' Function to calculate vector of max entropy noise distribution for a given value of b and d
#'
#' @param b 
#' @param d 
#' @param n 
#' @param tolentr 
#' @param tolq 
#' @param maxit 
#'
#' @return
#' @importFrom stats dpois
#' @author Johan Heldal
#'
#' @examples
calcprob <- function(b,
                     d,
                     n,
                     tolentr = 1.0e-6,
                     tolq = 1.0e-6,
                     maxit = 20) {
  if (d != b - 1 || n != b) {
    if (n <= d) {
      d1 <- -n
      d2 <- b - n
      x <- as.matrix(c(d1, d2:d))
    }
    if (d + 1 <= n) {
      d1 <- max(-d, b - n)
      x <- as.matrix(c(d1:d))
    }
    ##
    ## Calculate start values
    ##
    m <- nrow(x)
    en <- as.matrix(rep(1, length(x)))
    q0 <- sqrt(dpois(x - min(x), -min(x)))
    q0[1] <- -t(x[2:m, 1]) %*% q0[2:m] / x[1, 1]
    q0 <- q0 / sum(q0)
    E0x <- as.vector(t(x) %*% q0)
    sumq0 <- sum(q0)
    entropy0 <- sum(q0 * log(q0))
    l10 <- 1 + entropy0
    ql0 <- rbind(q0, l10)
    lagrange <- entropy0 + l10 * (sumq0 - 1)
    if (n != b + d) {
      l20 <- (sum(log(q0)) - entropy0) / sum(x)
      ql0 <- rbind(ql0, l20)
      lagrange <- lagrange + l20 * (E0x - 0)
    }
    ##
    ## Start iterations
    ##
    it <- 0
    dentr <- 100
    dq <- 100
    while ((dq > tolq) && (it <- it + 1) <= maxit) {
      ##
      ##
      Dlagr_Dq <- - (en + log(q0)) + l10 * en
      rl1 <- sum(q0) - 1
      r <- rbind(Dlagr_Dq, rl1)
      Jf1 <- cbind(-diag(as.vector(1 / q0)), en)
      Jf2 <- cbind(t(en), 0)
      Jf <- rbind(Jf1, Jf2)
      if (b + d != n) {
        Dlagr_Dq <- Dlagr_Dq + l20 * x
        rl2 <- t(x) %*% q0
        r <- rbind(r, rl2)
        Jf <- cbind(Jf, rbind(x, 0))
        Jf <- rbind(Jf, cbind(t(x), 0, 0))
      }
      dql <- solve(Jf, -r)
      ql1 <- ql0 + dql
      q1 <- ql1[1:m]
      q0 <- q1
      l10 <- ql1[m + 1, 1]
      l20 <- 0
      if (b + d != n)
        l20 <- ql0[m + 2, 1]
      ql0 <- ql1
      dq <- max(abs(as.vector(dql[1:m, 1])))
      entropy1 <- sum(q1 * log(q1))
      dentr < abs(entropy1 - entropy0)
      entropy0 <- entropy1
    }
  }
  if (d == b - 1 && n == b) {
    x <- as.matrix(0)
    q1 <- as.matrix(1)
    entropy1 <- 0
    l10 <- 0
    l20 <- 1
    it <- 0
    dq <- 0
    dentr <- 0
  }
  x2 <- x * x
  Varx <- t(x2) %*% q1
  Ex <- t(x) %*% q1
  Sump <- sum(q1)
  return(
    list(
      p = q1,
      X = x,
      entropy = -entropy1,
      lambda1 = l10,
      lambda2 = l20,
      iter = it,
      Sump = Sump,
      Ex = Ex,
      Varx = Varx,
      dentr = dentr,
      dq = dq
    )
  )
}
