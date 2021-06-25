#' Extract ORs for logistic regression models
#'
#' This function takes a regression object as an input and returns the
#' extracted values of beta coefficient, odds ratio (OR), 95% confidence
#' intervals (95% CIs), and P value for one or more given variable(s) of each model.
#'
#' @param regress_obj An object given by glm() function
#' @param index_lst A list of index(s) of variable(s) of interest
#' @param mark A mark created by users to identify model(s)
#' @return A table of extracted model values
#' @examples
#' # Create an example logistic model
#' example_model <- glm(y ~ x + c1 + c2, data = example_data, family = 'binomial')
#'
#' # Extract OR values for x
#' tab_1 <- extract_OR(example_model, 2, 'first_example')
#'
#' # Extract OR values for x and c1
#' tab_2 <- extract_OR(example_model, c(2, 3), 'second_example')
#' @export
extract_OR = function(regress_obj, index_lst, mark) {

  summary_reg = summary(regress_obj)
  coefficient = summary_reg$coefficients

  or_table = array(NA, dim = c(length(index_lst), 8))

  for(i in 1:length(index_lst)) {

    index = index_lst[i]
    OR = exp(coefficient[index])
    coeff = coefficient[index]

    confinterval_exp = exp(confint(regress_obj, level = 0.95))
    num_variable = length(confinterval_exp)/2
    lower_CI_exp = confinterval_exp[index]
    upper_CI_exp = confinterval_exp[num_variable + index]
    CI_95_exp = paste0('(', sprintf('%.2f', lower_CI_exp), ', ', sprintf('%.2f', upper_CI_exp), ')')

    confinterval = confint(regress_obj, level = 0.95)
    lower_CI = confinterval[index]
    upper_CI = confinterval[num_variable + index]
    CI_95 = paste0('(', sprintf('%.2f', lower_CI), ', ', sprintf('%.2f', upper_CI), ')')

    p = round(coefficient[num_variable*3 + index], 4)

    or_table[i,1] = mark
    or_table[i,2] = index
    or_table[i,3] = length(regress_obj$y)
    or_table[i,4] = sprintf('%.3f',coeff)
    or_table[i,5] = CI_95
    or_table[i,6] = round(OR, 2)
    or_table[i,7] = CI_95_exp
    or_table[i,8] = sprintf('%.3f',p)

    i = i + 1
  }

  or_table_final = as.data.frame(or_table)
  colnames(or_table_final) = c('mark','index', 'N', 'coeff', '95% CI', 'OR', '95% CI_OR', 'p-value')

  return(or_table_final)
}



#' Make an OR table for multiple logistic regression models
#'
#' This function takes a list of logistic regression models as an input and returns the
#' extracted values of beta coefficient, odds ratio (OR), 95% confidence
#' intervals (95% CIs), and P value for one or more given variable(s) of each model.
#'
#' @param lst_model A list of logistic regressionmodels
#' @param index_lst A list of index(s) of variable(s) of interest
#' @param model_head The mark created by users to identify model(s)
#' @return A table of extracted model values for multiple models
#' @examples
#' # Create an example logistic model lists
#' lst_model <- c('model_1', 'model_2', 'model_3')
#'
#' model_1 <- glm(y ~ x, data = example_data, family = 'binomial')
#' model_2 <- glm(y ~ x + c1, data = example_data, family = 'binomial')
#' model_3 <- glm(y ~ x + c1 + c2, data = example_data, family = 'binomial')
#'
#' # Make an OR table for variable x
#' eg_OR_tab <- make_OR_table(lst_model, 2, 'example_')
#' @export
make_OR_table = function(lst_model, index_lst, model_head) {

  # Make to a full table
  or_table = data.frame(
    Mark = character(),
    Index = integer(),
    N = integer(),
    coeff = double(),
    OR = double(),
    CI_95 = character(),
    p_value = character(),
    lower_CI = double(),
    upper_CI = double(),
    stringsAsFactors = F
  )

  for (i in 1:length(lst_model)) {

    model_object = get(lst_model[i])
    result = extract_OR(model_object, index_lst, paste0(model_head, lst_model[i]))

    or_table = rbind(or_table, result)

    i = i + 1
  }

  return(or_table)
}
