# Helper function: Validate and determine basis type
.validate_basis_type <- function(basis, model, basis_name = NULL) {
  type <- if (any(class(basis) %in% "crossbasis")) "cb" 
          else if (any(class(basis) %in% "onebasis")) "one" 
          else "gam"
  
  errormes <- "arguments 'basis' and 'model' not consistent. See help(crosspred)"
  
  if (type == "gam") {
    if (!is.character(basis) || length(basis) > 1L) stop(errormes)
    if (is.null(model) || !any(class(model) %in% "gam")) stop(errormes)
    
    name <- basis
    sterms <- sapply(model$smooth, function(x) x$term[1])
    
    if (name %in% sterms) {
      basis <- model$smooth[[which(sterms == name)[1]]]
    } else {
      stop(errormes)
    }
    
    if (length(which(sterms == name)) > 1) {
      warning(paste(name, "included in multiple smoothers, only the first one taken"))
    }
    
    if (!"cb.smooth" %in% class(basis) && basis$dim > 1L) {
      stop("predictions not provided for multi-dimensional smoothers other than 'cb'")
    }
  } else {
    name <- if (is.null(basis_name)) "basis" else basis_name
  }
  list(type = type, basis = basis, name = name)
}

# Helper function: Validate lag parameters
.validate_lag_params <- function(type, basis, lag, bylag, cumul) {
  origlag <- switch(type,
    cb = attr(basis, "lag"),
    one = c(0, 0),
    gam = if (is.null(basis$lag)) c(0, 0) else basis$lag
  )
  
  lag <- if (missing(lag)) origlag else mklag(lag)
  
  if (!all(lag == origlag) && cumul) {
    stop("cumulative prediction not allowed for lag sub-period")
  }
  
  lagfun <- switch(type,
    cb = attr(basis, "arglag")$fun,
    one = NULL,
    gam = if (basis$dim == 1L) NULL else basis$margin[[2]]$fun
  )
  
  if (bylag != 1L && !is.null(lagfun) && lagfun == "integer") {
    stop("prediction for non-integer lags not allowed for type 'integer'")
  }
  
  list(lag = lag, origlag = origlag)
}

# Helper function: Extract coefficients and variance-covariance matrix
.extract_coef_vcov <- function(type, basis, name, model, coef, vcov, model.link) {
  # Determine condition for coefficient selection
  cond <- switch(type,
    gam = with(basis, first.para:last.para),
    {
      if (ncol(basis) == 1L) {
        name
      } else if (type == "one") {
        paste0(name, "[[:print:]]*b[0-9]{1,2}")
      } else {
        paste0(name, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
      }
    }
  )
  
  if (!is.null(model)) {
    model.class <- class(model)
    coef <- getcoef(model, model.class)
    vcov <- getvcov(model, model.class)
    
    if (type == "gam") {
      indcoef <- indvcov <- cond
    } else {
      indcoef <- grep(cond, names(coef))
      indvcov <- grep(cond, rownames(vcov))
    }
    
    coef <- coef[indcoef]
    vcov <- vcov[indvcov, indvcov, drop = FALSE]
    model.link <- getlink(model, model.class, model.link)
  } else {
    model.class <- NA
  }
  
  # Validation
  npar <- if (type == "gam") length(indcoef) else ncol(basis)
  if (length(coef) != npar || length(coef) != nrow(vcov) || 
      any(is.na(coef)) || any(is.na(vcov))) {
    stop("coef/vcov not consistent with basis matrix. See help(crosspred)")
  }
  
  list(coef = coef, vcov = vcov, model.class = model.class, model.link = model.link)
}

# Helper function: Setup prediction variables and centering
.setup_prediction_vars <- function(type, basis, model, at, from, to, by, lag, bylag, cen) {
  range <- if (type == "gam") range(model$model[[basis$term[1]]]) else attr(basis, "range")
  
  at <- mkat(at, from, to, by, range, lag, bylag)
  predvar <- if (is.matrix(at)) rownames(at) else at
  predlag <- seqlag(lag, bylag)
  
  cen <- mkcen(cen, type, basis, range)
  
  # Clean up basis attributes
  if (type %in% c("one", "cb")) {
    if (type == "one") attributes(basis)$cen <- NULL
    if (type == "cb") attributes(basis)$argvar$cen <- NULL
  }
  
  list(at = at, predvar = predvar, predlag = predlag, cen = cen, range = range)
}


# Helper function: Calculate lag-specific effects
.calc_lag_effects <- function(type, basis, at, predvar, predlag, cen, coef, vcov) {
  Xpred <- mkXpred(type, basis, at, predvar, predlag, cen)
  
  matfit <- matrix(Xpred %*% coef, length(predvar), length(predlag))
  matse <- matrix(
    sqrt(pmax(0, rowSums((Xpred %*% vcov) * Xpred))), 
    length(predvar), length(predlag)
  )
  
  rownames(matfit) <- rownames(matse) <- predvar
  colnames(matfit) <- colnames(matse) <- outer("lag", predlag, paste, sep = "")
  
  list(matfit = matfit, matse = matse)
}


# Helper function: Calculate overall and cumulative effects
.calc_overall_effects <- function(type, basis, at, predvar, lag, cen, coef, vcov, cumul) {
  predlag <- seqlag(lag)
  Xpred <- mkXpred(type, basis, at, predvar, predlag, cen)
  
  Xpredall <- 0
  cumfit <- cumse <- NULL
  
  if (cumul) {
    cumfit <- cumse <- matrix(0, length(predvar), length(predlag))
  }
  
  for (i in seq(length(predlag))) {
    ind <- seq(length(predvar)) + length(predvar) * (i - 1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
    if (cumul) {
      cumfit[, i] <- Xpredall %*% coef
      cumse[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% vcov) * Xpredall)))
    }
  }
  
  allfit <- as.vector(Xpredall %*% coef)
  allse <- sqrt(pmax(0, rowSums((Xpredall %*% vcov) * Xpredall)))
  
  names(allfit) <- names(allse) <- predvar
  if (cumul) {
    rownames(cumfit) <- rownames(cumse) <- predvar
    colnames(cumfit) <- colnames(cumse) <- outer("lag", seqlag(lag), paste, sep = "")
  }
  
  list(allfit = allfit, allse = allse, cumfit = cumfit, cumse = cumse)
}

# Helper function: Add confidence intervals
.add_confidence_intervals <- function(results, model.link, ci.level, cumul) {
  z <- qnorm(1 - (1 - ci.level) / 2)
  is_exp_link <- !is.null(model.link) && model.link %in% c("log", "logit")
  
  # Matrix effects
  if (is_exp_link) {
    results$matRRfit <- exp(results$matfit)
    results$matRRlow <- exp(results$matfit - z * results$matse)
    results$matRRhigh <- exp(results$matfit + z * results$matse)
  } else {
    results$matlow <- results$matfit - z * results$matse
    results$mathigh <- results$matfit + z * results$matse
  }
  
  # Overall effects  
  if (is_exp_link) {
    results$allRRfit <- exp(results$allfit)
    results$allRRlow <- exp(results$allfit - z * results$allse)
    results$allRRhigh <- exp(results$allfit + z * results$allse)
    names(results$allRRlow) <- names(results$allRRhigh) <- names(results$allfit)
  } else {
    results$alllow <- results$allfit - z * results$allse
    results$allhigh <- results$allfit + z * results$allse
    names(results$alllow) <- names(results$allhigh) <- names(results$allfit)
  }
  
  # Cumulative effects
  if (cumul && !is.null(results$cumfit)) {
    if (is_exp_link) {
      results$cumRRfit <- exp(results$cumfit)
      results$cumRRlow <- exp(results$cumfit - z * results$cumse)
      results$cumRRhigh <- exp(results$cumfit + z * results$cumse)
    } else {
      results$cumlow <- results$cumfit - z * results$cumse
      results$cumhigh <- results$cumfit + z * results$cumse
    }
  }
  
  results
}

# Main function: Streamlined and modular
crosspred <- function(
    basis, model = NULL, coef = NULL, vcov = NULL, model.link = NULL, at = NULL,
    from = NULL, to = NULL, by = NULL, lag, bylag = 1, cen = NULL, ci.level = 0.95,
    cumul = FALSE) {
  
  # Input validation
  if (is.null(model) && (is.null(coef) || is.null(vcov))) {
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  }
  if (!is.numeric(ci.level) || ci.level >= 1 || ci.level <= 0) {
    stop("'ci.level' must be numeric and between 0 and 1")
  }
  
  # Process inputs
  basis_name <- deparse(substitute(basis))
  basis_info <- .validate_basis_type(basis, model, basis_name)
  lag_info <- .validate_lag_params(basis_info$type, basis_info$basis, lag, bylag, cumul)
  coef_info <- .extract_coef_vcov(basis_info$type, basis_info$basis, basis_info$name, 
                                  model, coef, vcov, model.link)
  pred_vars <- .setup_prediction_vars(basis_info$type, basis_info$basis, model, 
                                      at, from, to, by, lag_info$lag, bylag, cen)
  
  # Calculate effects
  lag_effects <- .calc_lag_effects(basis_info$type, basis_info$basis, pred_vars$at, 
                                   pred_vars$predvar, pred_vars$predlag, pred_vars$cen, 
                                   coef_info$coef, coef_info$vcov)
  overall_effects <- .calc_overall_effects(basis_info$type, basis_info$basis, pred_vars$at, 
                                           pred_vars$predvar, lag_info$lag, pred_vars$cen, 
                                           coef_info$coef, coef_info$vcov, cumul)
  
  # Assemble and finalize results
  results <- c(
    list(predvar = pred_vars$predvar, lag = lag_info$lag, bylag = bylag,
         coefficients = coef_info$coef, vcov = coef_info$vcov),
    lag_effects, overall_effects[c("allfit", "allse")],
    if (!is.null(pred_vars$cen)) list(cen = pred_vars$cen),
    if (cumul && !is.null(overall_effects$cumfit)) overall_effects[c("cumfit", "cumse")]
  )
  
  results <- .add_confidence_intervals(results, coef_info$model.link, ci.level, cumul)
  results <- c(results, list(ci.level = ci.level, model.class = coef_info$model.class, 
                            model.link = coef_info$model.link))
  
  class(results) <- "crosspred"
  results
}
