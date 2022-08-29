#' NHANES data from 2015-2016 and 2017-2018 cycles
#'
#' Combined NHANES data from the 2015-2016 and 2017-2018 cycles
#' The weights have been adjusted according to \url{https://wwwn.cdc.gov/nchs/nhanes/tutorials/module3.aspx}
#'
#' @format Data with the following variables:
#' \describe{
#'   \item{SEQN}{Respondent sequence number}
#'   \item{WTINT4YR}{Full sample 4 year interview weight}
#'   \item{WTMEC4YR}{Full sample 4 year MEC exam weight}
#'   \item{WTSB4YR}{Environmental B 4-year weights}
#'   \item{RIDSTATR}{Interview/Examination status}
#'   \item{RIAGENDR}{Gender of the participant}
#'   \item{RIDAGEYR}{Age in years of the participant at the time of screening.
#'   Individuals 80 and over are top-coded at 80 years of age.}
#'   \item{INDFMPIR}{A ratio of family income to poverty guidelines}
#'   \item{RIDRETH1}{Race/Hispanic origin}
#'   \item{INDHHIN2}{Total household income (reported as a range value in dollars)}
#'   \item{BMXBMI}{Body Mass Index (kg/m**2)}
#'   \item{BMXWAIST}{Waist Circumference (cm)}
#'   \item{BMXWT}{Weight (kg)}
#'   \item{BMXHT}{Standing Height (cm)}
#'   \item{URXUCR}{Creatinine, urine (mg/dL)}
#'   \item{URXCNP, URXCOP, URXECP, URXHIBP, URXMBP, URXMC1, URXMCOH, URXMEP,
#'   URXMHBP, URXMHH, URXMHNC, URXMHP, URXMIB, URXMNP, URXMOH, URXMZP}{Phthalates concentrations}
#'   \item{URDCNPLC, URDCOPLC, URDECPLC, URDHIBLC, URDMBPLC, URDMC1LC, URDMCOLC, URDMEPLC,
#'   URDMHBLC, URDMHHLC, URDMCHLC, URDMHPLC, URDMIBLC, URDMNPLC, URDMOHLC, URDMZPLC}{Phthalates comment code for whether the measurement is under the limit of detection}
#' }
#' @source \url{https://wwwn.cdc.gov/nchs/nhanes/search/default.aspx}
#' Detailed documentation of the phthalates variables can be found here:
#'   * \url{https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/PHTHTE_I.htm}
#'   * \url{https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/PHTHTE_J.htm}
"nhanes1518"
