lin <- function(x, intercept = FALSE) {
  nx <- names(x)
  x <- as.vector(x)
  # TRANSFORMATION
  basis <- as.matrix(x)
  if (intercept) basis <- cbind(1, basis)
  # NAMES AND ATTRIBUTES
  dimnames(basis) <- list(nx, seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis), list(intercept = intercept))
  class(basis) <- c("lin", "matrix")
  return(basis)
}

# scale用于标准化处理
poly <- function(x, degree = 1, scale, intercept = FALSE) {
  nx <- names(x)
  x <- as.vector(x)
  # TRANSFORMATION
  if (missing(scale)) scale <- max(abs(x), na.rm = TRUE)
  basis <- outer(x / scale, (1 - intercept):(degree), "^")

  # NAMES AND ATTRIBUTES
  dimnames(basis) <- list(nx, seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis), list(
    degree = degree, scale = scale,
    intercept = intercept
  ))
  class(basis) <- c("poly", "matrix")
  return(basis)
}

integer <- function(x, values, intercept = FALSE) {
  nx <- names(x)
  x <- as.vector(x)
  # DEFINE LEVELS AND TRANSFORM INTO A FACTOR
  levels <- if (!missing(values)) values else sort(unique(x))
  xfac <- factor(x, levels = levels)
  
  # TRANSFORMATION
  basis <- as.matrix(outer(xfac, levels, "==") + 0L)
  # IF INTERCEPT IS NOT REQUIRED, DROP THE FIRST COLUMN
  if (ncol(basis) > 1L) {
    if (!intercept) basis <- basis[, -1L, drop = FALSE] # 为何要去掉第一列？防止共线性
  } else {
    intercept <- TRUE
  }
  # NAMES AND ATTRIBUTES
  dimnames(basis) <- list(nx, seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis), list(
    values = levels,
    intercept = intercept
  ))
  class(basis) <- c("integer", "matrix")
  return(basis)
}
