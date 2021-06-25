#' Test for Mendelian randomization assumptions
#'
#' This function could be used for testing MR assumptions. Rules suggested by Glymour et al.
#' See details at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3366596/pdf/kwr323.pdf.
#'
#' @param gene_instrument The gene instrument used in the MR analysis (should be single SNP/polygenic score).
#' @param exposure Exposure of the MR model.
#' @param outcome Outcome of the MR model.
#' @param covariates Covariates adjusted in the model (default as NULL).
#' @param data_name Name of the dataset. Please notice that the dataset should have NO missing value.
#' @return A set of TRUE/FALSE indicating if the assumptions hold. If the MR assumptions hold,
#' all value should be TRUE.
#' @examples
#' MR_assumption_test('gene_instrument', 'exposure', 'outcome', covariates = 'a + b', 'example_dataset')
#' @export
MR_assumption_test <- function(gene_instrument, exposure, outcome, covariates = NULL, data_name) {

  # Create a set of regression
  if (!is.null(covariates)) {
    f1 = as.formula(paste(exposure, '~', gene_instrument, '+', covariates))
    f2 = as.formula(paste(outcome, '~', gene_instrument, '+', covariates))
    f3 = as.formula(paste(outcome, '~', exposure, '+', covariates))
    f4 = as.formula(paste(outcome, '~', exposure, '+', gene_instrument, '+', covariates))
  } else if (is.null(covariates)) {
    f1 = as.formula(paste(exposure, '~', gene_instrument))
    f2 = as.formula(paste(outcome, '~', gene_instrument))
    f3 = as.formula(paste(outcome, '~', exposure))
    f4 = as.formula(paste(outcome, '~', exposure, '+', gene_instrument))
  }

  f1_logreg = glm(f1, data = get(data_name), family = 'binomial')
  f2_logreg = glm(f2, data = get(data_name), family = 'binomial')
  f3_logreg = glm(f3, data = get(data_name), family = 'binomial')
  f4_logreg = glm(f4, data = get(data_name), family = 'binomial')

  summary_reg_1 = summary(f1_logreg)
  coefficient_1 = summary_reg_1$coefficients
  coef_gene_exposure = coefficient_1[2]

  summary_reg_2 = summary(f2_logreg)
  coefficient_2 = summary_reg_2$coefficients
  coef_gene_outcome = coefficient_2[2]

  summary_reg_3 = summary(f3_logreg)
  coefficient_3 = summary_reg_3$coefficients
  coef_exposure_outcome = coefficient_3[2]

  summary_reg_4 = summary(f4_logreg)
  coefficient_4 = summary_reg_4$coefficients
  coef_all_exposure = coefficient_4[2]
  coef_all_gene = coefficient_4[3]

  # Equation 1
  assumption_1 = coef_gene_outcome/coef_gene_exposure < coef_exposure_outcome
  print('Equation 1: The IV effect estimate is less than the ordinary least squares effect estimate')
  print(paste0(round(coef_gene_outcome/coef_gene_exposure, 4), ' < ', round(coef_exposure_outcome, 4)))
  print(assumption_1)
  print('==================================================================')
  # Equation 2
  assumption_2 = coef_all_gene < 0
  print('Equation 2: The coefficient of genetic instrument is negative')
  print(paste0(round(coef_all_gene, 4), ' < 0'))
  print(assumption_2)
  print('==================================================================')
  # Equation 3
  assumption_3 = abs(coef_exposure_outcome - coef_gene_outcome/coef_gene_exposure) < abs(coef_all_exposure - coef_gene_outcome/coef_gene_exposure)
  print('Equation 3: The causal slope adjusted for genetic instrument exceeds the bias of estimate without adjustment for genetic instrument')
  print(paste0('abs(', round(coef_exposure_outcome, 4), ' - ', round(coef_gene_outcome/coef_gene_exposure, 4), ') < abs(', round(coef_all_exposure, 4),' - ', round(coef_gene_outcome/coef_gene_exposure, 4), ')'))
  print(assumption_3)
  print('==================================================================')
  # Equation 4
  residual_MR = get(data_name)[outcome] - (coef_gene_outcome/coef_gene_exposure)*get(data_name)[exposure]
  cor_value = cor(residual_MR, get(data_name)[exposure])
  assumption_4 = cor_value > 0
  print('Equation 4: The residual is positively correlated with exposure')
  print(paste0(round(cor_value, 4), " > 0"))
  print(assumption_4[1])

}



#' Make a table to extract results from Mendelian randomization
#'
#' This function takes two glm() models as input: 1) exposure ~ gene_instrument + covariates;
#' 2) outcome ~ gene_instrument + covairates; and returns the extracted values of beta coefficient,
#' 95% confidence intervals (95% CIs), and P value for MR analysis.
#'
#' @param reg_object1 A regression model given by glm() function: exposure ~ gene_instrument (+ covariates)
#' @param reg_object2 A regression model given by glm() function: outcome ~ gene_instrument (+ covariates)
#' @param mark Mark of the MR model.
#' @return A table of extracted results from Mendelian randomization
#' @examples
#' # Create two models for MR analysis
#' f0_exposure <- glm(exposure ~ gene_instrument, data = example_data, family = 'binomial')
#' f0_outcome <- glm(outcome ~ gene_instrument, data = example_data, family = 'binomial')
#' crude_mr <- get_MR_value(f0_exposure, f0_outcome, 'example_MR')
#' @export
get_MR_value <- function(reg_object1, reg_object2, mark) {

  summary_reg_1 = summary(reg_object1)
  summary_reg_2 = summary(reg_object2)
  coefficient_1 = summary_reg_1$coefficients
  coefficient_2 = summary_reg_2$coefficients

  num_variables_reg1 = length(coefficient_1)/4
  se_position_reg1 = num_variables_reg1 + 2

  num_variables_reg2 = length(coefficient_2)/4
  se_position_reg2 = num_variables_reg2 + 2

  bx = coefficient_1[2]
  bxse = coefficient_1[se_position_reg1]
  by = coefficient_2[2]
  byse = coefficient_2[se_position_reg2]

  MRdata = MendelianRandomization::mr_input(bx = bx, bxse = bxse, by = by, byse = byse)

  IVWObject = MendelianRandomization::mr_ivw(MRdata,
                     model = "default",
                     robust = FALSE,
                     penalized = FALSE,
                     correl = FALSE,
                     weights = "simple",
                     psi = 0,
                     distribution = "normal",
                     alpha = 0.05)

  p_value = IVWObject$Pvalue
  estimate = IVWObject@Estimate
  CI_95 = paste0('(', sprintf('%.2f',IVWObject@CILower), ', ', sprintf('%.2f',IVWObject@CIUpper), ')')
  OR = exp(IVWObject@Estimate)
  lower = exp(IVWObject@CILower)
  upper = exp(IVWObject@CIUpper)
  CI_95_exp = paste0('(', sprintf('%.2f',lower), ', ', sprintf('%.2f',upper), ')')

  # Make to a series
  mr_series = c(mark, sprintf('%.2f',estimate), CI_95, sprintf('%.2f',OR), CI_95_exp, sprintf('%.3f',p_value))
  # Mark as a dataframe
  mr_df = as.data.frame(t(mr_series))
  colnames(mr_df) = c('mark', 'estimate', '95% CI', 'OR', '95% CI_OR', 'p-value')

  return(mr_df)
}
