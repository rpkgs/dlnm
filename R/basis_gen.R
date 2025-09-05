onebasis <- function(x, fun = "ns", ...) {
  nx <- names(x)
  x <- as.vector(x)
  range <- range(x, na.rm = TRUE)
  args <- list(...)
  args$x <- x
  cen <- args$cen

  # CHECK THE CONTENT AND OPTIONALLY MODIFY OBJECTS THROUGH assign
  checkonebasis(fun, args, cen)
  ###########################################################################
  # TRANSFORMATION
  # CREATE THE BASIS
  basis <- do.call(fun, args)
  # FORCE TO BE A MATRIX (NOT WITH as.matrix AS IT DELETES ATTRIBUTES)
  if (is.null(dim(basis))) dim(basis) <- c(length(x), 1)
  attr <- attributes(basis)
  ##########################################################################
  # NAMES AND ATTRIBUTES (KEEP cen IF PROVIDED, TO BE USED LATER FOR CENTERING)
  attributes(basis) <- c(list(fun = fun), attr, list(range = range, cen = cen))
  dimnames(basis) <- list(nx, paste("b", seq(ncol(basis)), sep = ""))
  class(basis) <- c("onebasis", "matrix")
  return(basis)
}


checkonebasis <- function(fun, args, cen) {
  # CHECK fun, AND IF fun HAS x ARGUMENT
  if (!is.character(fun)) stop("'fun' must be a string referring to a function")
  if (all(names(formals(fun)) != "x")) stop("'fun' must contain argument 'x'")
  # CHECK CENTERING, MOVED TO PREDICTION NOW, AND REMOVE FROM ARGUMENTS
  # ALSO, SET TO NULL FOR NON-CONTINUOUS FUNCTIONS
  if (!is.null(args$cen) && !"cen" %in% names(formals(fun))) {
    warning("centering through 'cen' now applied at the prediction stage. See ?crosspred")
    args$cen <- NULL
  }
  
  # OLD ARGUMENT bound FOR SPLINE FUNCTIONS
  if (fun %in% c("ns", "bs") && !is.null(args$bound)) {
    names(args)[names(args) == "bound"] <- "Boundary.knots"
    warning("use the default argument 'Boundary.knots' for fun 'ns'-'bs'")
  }

  # OLD THRESHOLD FUNCTIONS
  if (fun %in% c("hthr", "lthr", "dthr")) {
    args$side <- switch(fun,
      hthr = "h",
      lthr = "l",
      dthr = "d"
    )
    fun <- "thr"
    warning("function 'hthr'-'lthr'-'dthr' replaced by 'thr'. See ?thr")
  }
  if (fun == "thr" && !is.null(args$knots)) {
    names(args)[names(args) == "knots"] <- "thr.value"
    warning("argument 'knots' replaced by 'thr.value' in function thr. See ?thr")
  }
  # OLD STRATA FUNCTION
  if (fun == "strata" && !is.null(args$knots)) {
    names(args)[names(args) == "knots"] <- "breaks"
    warning("argument 'knots' replaced by 'breaks' in function strata. See ?strata")
  }
  assign("fun", fun, parent.frame())
  assign("args", args, parent.frame())
  assign("cen", cen, parent.frame())
}
