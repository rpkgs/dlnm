crossbasis <- function(x, lag, argvar = list(), arglag = list(), group = NULL, ...) {
  checkcrossbasis(argvar, arglag, list(...))

  lag <- if (missing(lag)) c(0, NCOL(x) - 1) else mklag(lag)

  ############################################################################
  # CREATE THE BASIS FOR THE PREDICTOR SPACE
  # x MUST BE A VECTOR OR MATRIX WITH NUMBER OF COLUMNS COMPATIBLE WITH lag
  # IF A VECTOR, x  IS TREATED AS A TIME SERIES
  # OTHERWISE, x IS TREATED AS A MATRIX OF LAGGED OCCURRENCES
  x <- as.matrix(x)
  dim <- dim(x)
  if (!dim[2] %in% c(1L, diff(lag) + 1L)) {
    stop(
      "NCOL(x) must be equal to 1 (if x is a time series vector), ",
      "otherwise to the lag period (for x as a matrix of lagged occurrences)"
    )
  }

  # THE BASIS TRANSFORMATION CREATES DIFFERENT MATRICES DEPENDING THE DATA :
  #   IF TIME SERIES, EACH COLUMN CONTAINS THE UNLAGGED TRANSFORMATION
  #   IF NOT, EACH COLUMN CONTAINS THE TRANFORMATION FOR ALL THE LAGGED VALUES
  basisvar <- do.call("onebasis", modifyList(argvar, list(x = as.numeric(x))))

  ############################################################################
  # CREATE THE BASIS FOR THE LAG SPACE
  # SET FUN="STRATA" AND DF=1 UNDER SPECIFIC CIRCUMSTANCES
  if (length(arglag) == 0L || diff(lag) == 0L) {
    arglag <- list(fun = "strata", df = 1, intercept = TRUE)
  }

  # IF NOT SPECIFIED AND AN ARGUMENT, INCLUDE AN INTERCEPT BY DEFAULT
  if ((is.null(arglag$fun) || "intercept" %in% names(formals(arglag$fun))) &&
    sum(pmatch(names(arglag), "intercept", nomatch = 0)) == 0) {
    arglag$intercept <- TRUE
  }

  # FORCE UNCENTERED TRANSFORMATIONS
  arglag$cen <- NULL
  # THE BASIS TRANSFORMATIONS ARE ONLY APPLIED TO THE LAG VECTOR
  # DIMENSIONS ACCOUNTED FOR IN CROSS-BASIS COMPUTATIONS BELOW
  basislag <- do.call("onebasis", modifyList(arglag, list(x = seqlag(lag))))

  ############################################################################
  # CROSSBASIS COMPUTATION
  # GROUP
  if (!is.null(group)) checkgroup(group, x, basisvar, lag)
  # COMPUTE CROSS-BASIS:
  #   FOR TIME SERIES DATA, COMPUTE THE MATRIX OF LAGGED OCCURRENCES FIRST
  #   IF x WAS ALREADY A MATRIX, JUST RECOMPUTE THE APPROPRIATE DIMENSIONS
  #   NB: ORDER OF TRANSFORMATION IN THE TENSOR CHANGED SINCE VERSION 2.2.4
  vx <- ncol(basisvar) # var: [n, vx]
  vl <- ncol(basislag) # lag: [nlag, vl]

  crossbasis <- matrix(0, nrow = dim[1], ncol = vx * vl)
  for (v in seq(length = vx)) { # for t
    if (dim[2] == 1L) {
      mat <- as.matrix(Lag(basisvar[, v], seqlag(lag), group = group))
      # mat: {X[n, t], X[n, t-1], ..., X[n, t-L]}
    } else {
      mat <- matrix(basisvar[, v], ncol = diff(lag) + 1)
    }
    for (l in seq(length = vl)) {
      ck <- basislag[, l] #
      crossbasis[, vl * (v - 1) + l] <- mat %*% ck # 精彩
    }
  }

  ############################################################################
  # ATTRIBUTES AND NAMES
  # NAMES
  #   NB: ORDER CHANGED SINCE VERSION 2.2.4
  cn <- paste0("v", rep(seq(vx), each = vl), ".l", rep(seq(vl), vx))
  dimnames(crossbasis) <- list(rownames(x), cn)
  # REDEFINE ARGUMENTS FOR BASES, THEY MIGHT HAVE BEEN CHANGED BY onebasis FIRST VAR
  ind <- match(names(formals(attributes(basisvar)$fun)),
    names(attributes(basisvar)),
    nomatch = 0
  )
  argvar <- c(attributes(basisvar)["fun"], attributes(basisvar)[ind])
  # THEN LAG
  ind <- match(names(formals(attributes(basislag)$fun)),
    names(attributes(basislag)),
    nomatch = 0
  )
  arglag <- c(attributes(basislag)["fun"], attributes(basislag)[ind])
  # THEN ADD CENTERING FOR VAR, IF PROVIDED (OTHERWISE NULL)
  argvar$cen <- attributes(basisvar)$cen

  attributes(crossbasis) <- c(attributes(crossbasis), list(
    df = c(vx, vl), range = range(x, na.rm = T), lag = lag,
    argvar = argvar, arglag = arglag, 
    basisvar = basisvar, basislag = basislag
  ))
  if (!is.null(group)) attributes(crossbasis)$group <- length(unique(group))

  class(crossbasis) <- c("crossbasis", "matrix")
  return(crossbasis)
}


checkgroup <- function(group, x, basisvar, lag) {
  if (NCOL(x) > 1L) stop("'group' allowed only for time series data")
  if (min(tapply(x, group, length)) <= diff(lag)) {
    stop("each group must have length > diff(lag)")
  }
}

checkcrossbasis <- function(argvar, arglag, addarg) {
  # CHECK LIST FORMAT
  if (!is.list(argvar)) stop("'argvar' must be a list")
  if (!is.list(arglag)) stop("'arglag' must be a list")
  #
  # OLD ARGUMENT type
  if (is.null(argvar$fun) && !is.null(argvar$type)) {
    names(argvar)[names(argvar) == "type"] <- "fun"
    assign("argvar", argvar, parent.frame())
    warning("argument 'type' replaced by 'fun'. See ?onebasis")
  }
  if (is.null(arglag$fun) && !is.null(arglag$type)) {
    names(arglag)[names(arglag) == "type"] <- "fun"
    assign("arglag", arglag, parent.frame())
    warning("argument 'type' replaced by 'fun'. See ?onebasis")
  }
  # 'VERY' OLD USAGE
  if (any(c(
    "vartype", "vardf", "vardegree", "varknots", "varbound", "varint",
    "cen", "cenvalue", "maxlag", "lagtype", "lagdf", "lagdegree", "lagknots",
    "lagbound", "lagint"
  ) %in% addarg)) {
    stop("old usage not allowed any more. See ?crossbasis and ?onebasis")
  }
}
