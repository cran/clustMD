#' Model based clustering for mixed data: clustMD
#'
#' Model-based clustering of mixed data (i.e. data that consist of continuous,
#' binary, ordinal or nominal variables) using a parsimonious mixture of latent
#' Gaussian variable models.
#' @aliases clustMD-package
#' @author Damien McParland
#' 
#' Damien McParland <damien.mcparland@ucd.ie>
#' Isobel Claire Gormley <claire.gormley@ucd.ie>
#' 
#' @references McParland, D. and Gormley, I.C. (2016). Model based clustering 
#'     for mixed data: clustMD. Advances in Data Analysis and Classification, 
#'     10 (2):155-169.
#'
#' @seealso \code{\link{clustMD}}
#'
#' @docType package
#' @keywords package
#' 
"_PACKAGE"

# ------------------------------------------------------------------- #

#################
### Data Sets ###
#################

#' Byar prostate cancer data set.
#'
#' A data set consisting of variables of mixed type measured on a group of
#' prostate cancer patients. Patients have either stage 3 or stage 4 prostate
#' cancer.
#' 
#' @format  A data frame with 475 observations on the following 15 variables.
#' \describe{
#'  \item{\code{Age}}{a numeric vector indicating the age of the patient.}
#'  \item{\code{Weight}}{a numeric vector indicating the weight of the patient.}
#'  \item{\code{Performance.rating}}{an ordinal variable indicating how active
#'      the patient is: 0 - normal activity, 1 - in bed less than 50\% of 
#'      daytime, 2 - in bed more than 50\% of daytime, 3 - confined to bed.}
#'  \item{\code{Cardiovascular.disease.history}}{a binary variable indicating
#'      if the patient has a history of cardiovascular disease: 0 - no, 1 - 
#'      yes.}
#'  \item{\code{Systolic.Blood.pressure}}{a numeric vector indicating the 
#'      systolic blood pressure of the patient in units of ten.}
#'  \item{\code{Diastolic.blood.pressure}}{a numeric vector indicating the 
#'      diastolic blood pressure of the patient in units of ten.}
#'  \item{\code{Electrocardiogram.code}}{a nominal variable indicating the 
#'      electorcardiogram code: 0 - normal, 1 - benign, 2 - rythmic 
#'      disturbances and electrolyte changes, 3 - heart blocks or conduction
#'      defects, 4 - heart strain, 5 - old myocardial infarct, 6 - recent 
#'      myocardial infarct.}
#'  \item{\code{Serum.haemoglobin}}{a numeric vector indicating the serum
#'      haemoglobin levels of the patient measured in g/100ml.}
#'  \item{\code{Size.of.primary.tumour}}{a numeric vector indicating the 
#'      estimated size of the patient's primary tumour in centimeters squared.}
#'  \item{\code{Index.of.tumour.stage.and.histolic.grade}}{a numeric vector 
#'      indicating the combined index of tumour stage and histolic grade of the
#'      patient.}
#'  \item{\code{Serum.prostatic.acid.phosphatase}}{a numeric vector indicating 
#'      the serum prostatic acid phosphatase levels of the patient in 
#'      King-Armstong units.}
#'  \item{\code{Bone.metastases}}{a binary vector indicating the presence of 
#'      bone metastasis: 0 - no, 1 - yes.}
#'  \item{\code{Stage}}{the stage of the patient's prostate cancer.}
#'  \item{\code{Observation}}{a patient ID number.}
#'  \item{\code{SurvStat}}{the post trial survival status of the patient: 
#'      0 - alive, 1 - dead from prostatic cancer, 2 - dead from heart or 
#'      vascular disease, 3 - dead from cerebrovascular accident, 3 - dead form
#'      pulmonary ebolus, 5 - dead from other cancer, 6 - dead from respiratory
#'      disease, 7 - dead from other specific non-cancer cause, 8 - dead from 
#'      other unspecified non-cancer cause, 9 - dead from unknown cause.}
#'}
#'
#' @source Byar, D.P. and Green, S.B. (1980). The choice of treatment for 
#'     cancer patients based on covariate information: applications to prostate
#'     cancer. Bulletin du Cancer 67: 477-490.
#'     
#'     Hunt, L., Jorgensen, M. (1999). Mixture model clustering using the 
#'     multimix program. Australia and New Zealand Journal of Statistics 41:
#'     153-171.
#'
#' @keywords datasets
#' 
"Byar"