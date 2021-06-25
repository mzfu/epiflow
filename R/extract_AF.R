#' Extract a table for population attributable fraction for logistic regression models
#'
#' This function takes a logistic regression object as an input and returns an
#' object of PAF values and its 95% confidence intervals (95% CIs).
#'
#' @param object A regression object given by glm() function
#' @param data An optional data frame, list or environment (or object coercible by as.data.frame
#' to a data frame) containing the variables in the model. If not found in data, the
#' variables are taken from environment (formula), typically the environment from
#' which the function is called.
#' @param exposure The name of the exposure variable as a string.
#' @param mark A mark created by the user to identify model(s)
#' @param case.control Can be set to TRUE if the data is from a non-matched case control study. By
#' default case.control is set to FALSE which is used for cross-sectional sampling designs.
#' @param confidence.level User-specified confidence level for the confidence intervals. If not specified it
#' defaults to 95 percent. Should be specified in decimals such as 0.95 for 95 percent.
#' @return A table of values for PAF.
#' @examples
#' # Create an example logistic model
#' example_model <- glm(y ~ x + c1 + c2, data = example_data, family = 'binomial')
#' return_AF <- extract_AF(object = example_model, data = example_data, exposure = "x", mark = 'example_1')
#' @export
extract_AF <- function(object, data, exposure, mark, case.control = FALSE, confidence.level = 0.95) {

  call = match.call()
  # Warning if the object is not a glm object
  if(!(as.character(object$call[1]) == "glm"))
    stop("The object is not a glm object", call. = FALSE)
  # Warning if the object is not a logistic regression
  if(!(object$family[1] == "binomial" & object$family[2] == "logit"))
    stop("The object is not a logistic regression", call. = FALSE)
  #### Preparation of dataset ####
  formula = object$formula
  #data <- object$data
  npar = length(object$coef)

  ## Delete rows with missing on variables in the model ##
  data = as.data.frame(data)
  rownames(data) = 1:nrow(data)
  m = model.matrix(object = formula, data = data)
  complete = as.numeric(rownames(m))
  data = data[complete, ]
  outcome = as.character(terms(formula)[[2]])
  n = nrow(data)
  n.cases = sum(data[, outcome])
  n.cluster <- 0


  ## Checks ##
  # if(!is.binary(data[, exposure]))
  #   stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)

  # Create dataset data0 for counterfactual X = 0
  data0 = data
  data0[, exposure] = 0

  ## Design matrices ##
  design = model.matrix(object = delete.response(terms(object)), data = data)
  design0 = model.matrix(object = delete.response(terms(object)), data = data0)

  #### Meat: score equations ####
  ## If sampling design is case-control ##
  if (case.control == TRUE){
    ## Create linear predictors to estimate the log odds ratio ##
    diff.design = design0 - design
    linearpredictor = design  %*% coef(object)
    linearpredictor0 = design0 %*% coef(object)
    #log odds ratio#
    log.or = linearpredictor - linearpredictor0
    ## Estimate approximate AF ##
    AF.est   = 1 - sum(data[, outcome] * exp( - log.or)) / sum(data[, outcome])
    #### Meat: score equations ####
    ## Score equation 1 ## individual estimating equations of the estimate of AF
    score.AF = data[, outcome] * (exp( - log.or) - AF.est)
    ## Score equation 2 ## individual estimating equations from conditional logistic reg.
    pred.diff = data[, outcome] - predict(object, newdata = data, type = "response")
    score.beta = design * pred.diff
    score.equations = cbind(score.AF, score.beta)
    meat <- var(score.equations, na.rm=TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    #### Estimating variance using Sandwich estimator ####
    hessian.AF1 = - data[, outcome]
    hessian.AF2 = (design0 - design) * as.vector(data[, outcome] * exp( - log.or))
    hessian.AF = cbind(mean(hessian.AF1), t(colMeans(hessian.AF2, na.rm = TRUE)))
    hessian.beta = cbind(matrix(rep(0, npar), nrow = npar, ncol = 1), - solve(vcov(object = object)) / n)
    ### Bread ###
    bread = rbind(hessian.AF, hessian.beta)
    #### Sandwich ####
    sandwich = (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    AF.var = sandwich[1, 1]
    #### Output ####
    out = c(list(AF.est = AF.est, AF.var = AF.var, log.or = log.or,
                 objectcall = object$call, call = call, exposure = exposure, outcome = outcome, object = object,
                 sandwich = sandwich, formula = formula,
                 n = n, n.cases = n.cases, n.cluster = n.cluster))
  }
  ## If sampling design is cross-sectional ##
  else {
    ## Score equation 1 ##
    score.P = data[, outcome]
    pred.Y  = predict(object, newdata = data, type = "response")
    ## Score equation 2 ##
    score.P0 = predict(object, newdata = data0, type = "response")
    ## Score equation 3 ##
    score.beta = design * (score.P - pred.Y)
    ### Meat ###
    score.equations = cbind(score.P, score.P0, score.beta)
    meat = var(score.equations, na.rm = TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    hessian.P = matrix(c(- 1, 0, rep(0,npar)), nrow = 1, ncol = 2 + npar)
    ## Hessian of score equation 2 ##
    g = family(object)$mu.eta
    dmu.deta = g(predict(object = object, newdata = data0))
    deta.dbeta = design0
    dmu.dbeta = dmu.deta * deta.dbeta
    hessian.P0 = matrix(c(0, - 1, colMeans(dmu.dbeta)), nrow = 1, ncol = 2 + npar)
    ## Hessian of score equation 3 ##
    hessian.beta = cbind(matrix(rep(0, npar * 2), nrow = npar, ncol = 2)
                         , - solve(vcov(object = object)) / n)
    ### Bread ###
    bread = rbind(hessian.P, hessian.P0, hessian.beta)
    #### Sandwich ####
    sandwich = (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    #### Point estimate of AF ####
    P.est  = mean(score.P, na.rm = TRUE)
    P0.est = mean(score.P0, na.rm = TRUE)
    AF.est = 1 - P0.est / P.est
    ## Delta method for variance estimate ##
    gradient = as.matrix(c(P0.est / P.est ^ 2, - 1 / P.est), nrow = 2, ncol = 1)
    AF.var = t(gradient) %*% sandwich %*% gradient
    P.var = sandwich[1, 1]
    P0.var = sandwich[2, 2]

    objectcall = object$call
    #### Output ####
    out = c(list(AF.est = AF.est, AF.var = AF.var, P.est = P.est, P0.est = P0.est, P.var = P.var,
                 P0.var = P0.var, objectcall = objectcall, call = call, exposure = exposure, outcome = outcome,
                 object = object, sandwich = sandwich, gradient = gradient, formula = formula,
                 n = n, n.cases = n.cases, n.cluster = n.cluster))
  }

  AF_table <- array(NA, dim=c(1,6))

  AF_est <- AF.est
  N <- n
  se <- sqrt(AF.var)
  zvalue <- AF.est / sqrt(AF.var)
  p <- 2 * pnorm( - abs(zvalue))
  lower <- AF.est - abs(qnorm((1 - confidence.level) / 2)) * se
  upper <- AF.est + abs(qnorm((1 - confidence.level) / 2)) * se

  AF_table[1,1] = mark
  AF_table[1,2] = N
  AF_table[1,3] = sprintf('%.4f', AF_est)
  AF_table[1,4] = sprintf('%.4f', lower)
  AF_table[1,5] = sprintf('%.4f', upper)
  AF_table[1,6] = sprintf('%.3f', p)

  AF_table_final = as.data.frame(AF_table)
  colnames(AF_table_final) = c('mark', 'N', 'AF', 'lower', 'upper', 'p-value')

  return(AF_table_final)
}
