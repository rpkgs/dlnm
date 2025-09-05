coef.crossreduce <- function(object, ...) {
  return(object$coef)
}

coef.crosspred <- function(object, ...) {
  return(object$coef)
}

vcov.crossreduce <- function(object, ...) {
  return(object$vcov)
}

vcov.crosspred <- function(object, ...) {
  return(object$vcov)
}


getcoef <- function(model, class) {
  # NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
  coef <- if (any(class %in% c("glm", "gam", "coxph"))) {
    coef(model)
  } else if (any(class %in% c("lme", "lmerMod", "glmerMod", "lmerModLmerTest"))) {
    fixef(model)
  } else {
    tryCatch(coef(model), error = function(w) "error")
  }
  if (identical(coef, "error")) {
    stop(
      "methods for coef() and vcov() must ",
      "exist for the class of object 'model'. If not, extract them manually and ",
      "use the arguments 'coef' and 'vcov'"
    )
  }
  return(coef)
}

getlink <- function(model, class, model.link = NULL) {
  # IF PROVIDED, JUST RETURN
  if (!is.null(model.link)) {
    return(model.link)
  }
  # OTHERWISE, EXTRACT FROM MODEL (IF AVAILABLE)
  link <- if (all(class %in% c("lm")) || all(class %in% c("lme")) ||
    any(class %in% "nlme") || any(class %in% "lmerMod")) {
    "identity"
  } else if (any(class %in% c("clogit"))) {
    "logit"
  } else if (any(class %in% c("coxph"))) {
    "log"
  } else if (any(class %in% c("glm")) || any(class %in% c("glmmPQL"))) {
    model$family$link
  } else if (any(class %in% c("glmerMod"))) {
    model@resp$family$link
  } else {
    NA
  }
  return(link)
}

getvcov <- function(model, class) {
  # EXTRACT VCOV
  # NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
  vcov <- if (any(class %in% c("lm", "glm", "lme", "coxph")) &&
    !identical(class, c("gee", "glm"))) {
    vcov(model)
  } else if (identical(class, c("gee", "glm"))) {
    model$robust.variance
  } else if (any(class %in% c("lmerMod", "glmerMod", "lmerModLmerTest"))) {
    as.matrix(vcov(model))
  } else {
    tryCatch(vcov(model), error = function(w) "error")
  }
  if (identical(vcov, "error")) {
    stop(
      "methods for coef() and vcov() must ",
      "exist for the class of object 'model'. If not, extract them manually and ",
      "use the arguments 'coef' and 'vcov'"
    )
  }
  return(vcov)
}
