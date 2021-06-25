#' Extract coefficients for linear regression models
#'
#' This function takes a regression object as an input and returns the
#' extracted values of beta coefficient, 95% confidence intervals (95% CIs),
#' and P value for one or more given variable(s) of each model.
#'
#' @param regress_obj An object given by lm() function
#' @param index_lst A list of index(s) of variable(s) of interest
#' @param mark A mark created by users to identify model(s)
#' @return A table of extracted model values
#' @examples
#' # Create an example logistic model
#' example_model <- lm(y ~ x + c1 + c2, data = example_data)
#'
#' # Extract OR values for x
#' tab_1 <- extract_coeff_lm(example_model, 2, 'first_example')
#'
#' # Extract OR values for x and c1
#' tab_2 <- extract_coeff_lm(example_model, c(2, 3), 'second_example')
#' @export
extract_coeff_lm = function(regress_obj, index_lst, mark) {

  summary_reg = summary(regress_obj)
  coefficient = summary_reg$coefficients

  coef_table = array(NA, dim=c(length(index_lst), 6))

  for(i in 1:length(index_lst)) {

    index = index_lst[i]
    coeff = coefficient[index]
    confinterval = confint(regress_obj, level = 0.95)
    num_variable = length(confinterval)/2
    lower_CI = confinterval[index]
    upper_CI = confinterval[num_variable + index]
    CI_95 = paste0('(', sprintf('%.2f',lower_CI), ', ', sprintf('%.2f',upper_CI), ')')

    p = round(coefficient[num_variable*3 + index], 4)

    coef_table[i,1] = mark
    coef_table[i,2] = index
    coef_table[i,3] = length(summary_reg$residuals)
    coef_table[i,4] = sprintf('%.3f',coeff)
    coef_table[i,5] = CI_95
    coef_table[i,6] = sprintf('%.3f',p)

    i = i + 1
  }

  coef_table_final = as.data.frame(coef_table)
  colnames(coef_table_final) = c('mark','index', 'N', 'coeff', '95% CI', 'p-value')

  return(coef_table_final)
}



#' Make a coefficient table for multiple linear regression models
#'
#' This function takes a list of linear regression models as an input and returns the
#' extracted values of beta coefficient, 95% confidence intervals (95% CIs),
#' and P value for one or more given variable(s) of each model.
#'
#' @param lst_model A list of linear regressionmodels
#' @param index_lst A list of index(s) of variable(s) of interest
#' @param model_head The mark created by users to identify model(s)
#' @return A table of extracted model values for multiple models
#' @examples
#' # Create an example linear model lists
#' lst_model <- c('model_1', 'model_2', 'model_3')
#'
#' model_1 <- lm(y ~ x, data = example_data)
#' model_2 <- lm(y ~ x + c1, data = example_data)
#' model_3 <- lm(y ~ x + c1 + c2, data = example_data)
#'
#' # Make a coefficient table for variable x
#' eg_linear_tab <- make_coeff_table_lm(lst_model, 2, 'example_')
#' @export
make_coeff_table_lm = function(lst_model, index_lst, model_head) {

  # Make to a full table
  coeff_table = data.frame(
    Mark = character(),
    Index = integer(),
    N = integer(),
    coeff = double(),
    CI = double(),
    p_value = double(),
    stringsAsFactors = F
  )

  for (i in 1:length(lst_model)) {

    model_object = get(lst_model[i])
    try = summary(model_object)

    result = extract_coeff_lm(model_object, index_lst, paste0(model_head, lst_model[i]))
    coeff_table = rbind(coeff_table, result)

    i = i + 1
  }

  return(coeff_table)
}
