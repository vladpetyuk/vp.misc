
#' SRM proteomics data from human CSF
#'
#' Targeted selected reaction monitoring (SRM) analysis of proteins in
#' human cerebrospinal fluid.
#'
#' @docType data
#'
#' @usage data(srm_msnset)
#'
#' @format An \code{MSnSet} with 42 features and 42 samples.
#'
#' @section phenoData:
#' A data frame consisting of 42 rows (samples) and 82 columns:
#' \describe{
#'   \item{sample.id}{A unique sample identifier. The match group and subject
#'   type.}
#'   \item{SampleNum}{A unique sample identifier. Same as \code{Prep Alias}.}
#'   \item{match.group}{Strings "01", "02", "04", "05", ..., "19".}
#'   \item{subject.type}{Either "case", "control.1", or "control.2".}
#'   \item{projid}{A unique sample identifier.}
#'   \item{age_death}{Age at death.}
#'   \item{msex}{Sex of the subject. Male is denoted by 1, female by 0.}
#'   \item{race}{The race of the subject. All entries are "White".}
#'   \item{spanish}{All entries are "Not Hispanic".}
#'   \item{educ}{}
#'   \item{pmi}{An estimate of the time since death (postmortem interval).}
#'   \item{parksc_lv}{}
#'   \item{dxpark}{}
#'   \item{cog_reduced}{}
#'   \item{dlbdx}{}
#'   \item{henl_4gp}{}
#'   \item{niareagansc}{}
#'   \item{amelidrg_ever}{}
#'   \item{lb_dens_lc:lb_dens}{}
#'   \item{lb_dens_sn_state}{Non-missing entries are "zero" or "high".}
#'   \item{tang_dens_sn}{}
#'   \item{tang_dens_sn_state}{Non-missing entries are "zero", "low", or "mid".}
#'   \item{amy_mean_sn}{}
#'   \item{amy_mean_sn_state}{All non-missing entries are "low".}
#'   \item{dlbany}{}
#'   \item{dlb_acc}{}
#'   \item{dlb_ag}{}
#'   \item{dlb_ec}{}
#'   \item{dlb_mfc}{}
#'   \item{dlb_mtc}{}
#'   \item{dlb_sn}{}
#'   \item{acaccb12}{}
#'   \item{acing12}{}
#'   \item{ah12}{}
#'   \item{aitb12}{}
#'   \item{amfc12}{}
#'   \item{amyag}{}
#'   \item{amyec}{}
#'   \item{asf12}{}
#'   \item{tauag}{}
#'   \item{tauec}{}
#'   \item{tcaccb12}{}
#'   \item{tcing12}{}
#'   \item{th12}{}
#'   \item{titb12}{}
#'   \item{tmfc12}{}
#'   \item{tsf12}{}
#'   \item{amy_mean_sn.predicted}{}
#'   \item{amy_mean_sn.imputed}{}
#'   \item{amy_mean_state.predicted}{All non-missing entries are "low".}
#'   \item{amy_mean_state.imputed}{All entries are "low".}
#'   \item{tang_dens_sn.predicted}{}
#'   \item{tang_dens_sn.imputed}{}
#'   \item{tang_dens_sn_state.predicted}{One of "zero", "low", "mid" or "high".}
#'   \item{tang_dens_sn_state.imputed}{One of "zero", "low", "mid" or "high".}
#'   \item{lb_dens_sn.predicted}{}
#'   \item{lb_dens_sn.imputed}{}
#'   \item{lb_dens_sn_state.predicted}{Either "zero" or "high".}
#'   \item{lb_dens_sn_state.imputed}{Either "zero" or "high".}
#'   \item{inclusion}{All entries are \code{TRUE}.}
#'   \item{case}{\code{TRUE} if \code{subject.type} is "case".}
#'   \item{control}{\code{TRUE} if \code{subject.type} is "control.1" or
#'   "control.2".}
#'   \item{match.1}{}
#'   \item{match.2}{}
#'   \item{slp_resid}{}
#'   \item{Prep order}{Integers 1-3.}
#'   \item{Prep Alias}{A unique sample identifier. Same as \code{SampleNum}.}
#'   \item{PrepBatch}{Batch number within \code{Prep order}. Integers 1-3.}
#'   \item{Coomassie (pre-digest)}{}
#'   \item{Visible blood contamination}{"contaminated" or blank ("").}
#'   \item{SPE Batch}{Integers 1-4.}
#'   \item{BCA (post-digest)}{The BCA Protein Assay concentration in mg/mL?
#'   Used for total protein quantitation.}
#'   \item{Yield (ug)}{}
#'   \item{Visible blood contamination after SPE}{"contaminated" or blank ("").}
#'   \item{Dilution Volume (Volume Remaing)}{}
#'   \item{dataset name}{The dataset name. The format is "Levy_CSF_SRM_XX",
#'   where "XX" is the \code{Prep Alias}.}
#' }
#'
#' @source Add source.
#'
#' @examples
#' data("srm_msnset")
#' head(pData(msnset))
#' head(fData(msnset))
#' head(exprs(msnset))
"msnset"


#' CPTAC Ovarian Cancer Proteomics Dataset
#'
#' Processed iTRAQ data from both PNNL and JHU sites.
#'
#' @docType data
#'
#' @usage data("cptac_oca")
#'
#' @format An \code{MSnSet} with 8103 features and 73 samples.
#'
#' @section pData:
#' A data frame consisting of 73 rows (samples) and 39 columns:
#' \describe{
#'   \item{Tumor.Nuclei}{Percent of cells in the tissue sample with tumor
#'   nuclei.}
#'   \item{Necrosis}{Percent of tumor nuclei that are necrotic.}
#'   \item{PNNL.Weight..mg.}{Weight of the non-tumor tissue in milligrams.}
#'   \item{Residual.Tumor.Tissue..mg.}{Weight of the tumor tissue in
#'   milligrams.}
#'   \item{TSS}{}
#'   \item{primary_therapy_outcome_success}{One of "COMPLETE.RESPONSE",
#'   "STABLE.DISEASE", or "PROGRESSIVE.DISEASE".}
#'   \item{days_to_death}{Days between }
#'   \item{Platinum.Refractory}{}
#'   \item{tumor_grade}{}
#'   \item{tumor_residual_disease}{}
#'   \item{tumor_stage}{Tumor stage.}
#'   \item{anatomic_organ_subdivision}{Either "left", "right", or "bilateral".}
#'   \item{days_to_birth}{Negative age (in days) of subject.}
#'   \item{days_to_tumor_recurrence}{Days before first tumor recurrence.}
#'   \item{race}{The race of the subject.}
#'   \item{site_of_tumor_first_recurrence}{"LOCO.REGIONAL" or "METASTASIS".}
#'   \item{Batch}{The sample batch.}
#'   \item{additional_chemo_therapy}{"YES" if the subject received additional
#'   chemotherapy; "NO" otherwise.}
#'   \item{additional_drug_therapy}{"YES" if the subject received additional
#'   drug therapy; "NO" otherwise.}
#'   \item{additional_pharmaceutical_therapy}{"YES" if the subject received
#'   additional pharmaceutical therapy; "NO" otherwise.}
#'   \item{additional_radiation_therapy}{"YES" if the subject received
#'   additional radiation therapy; "NO" otherwise.}
#'   \item{targeted_molecular_therapy}{"YES" if the subject received
#'   targeted molecular therapy; "NO" otherwise.}
#'   \item{year_of_initial_pathologic_diagnosis}{The year when subjects were
#'   first diagnosed with ovarian cancer.}
#'   \item{iTRAQ_Batch}{The iTRAQ batch identifier.}
#'   \item{iTRAQ_ID}{}
#'   \item{ReporterIon}{}
#'   \item{SUBTYPE}{Ovarian cancer subtype. Either "Differentiated",
#'   "Immunoreactive", "Mesenchymal", or "Proliferative".}
#'   \item{AGE}{Subject age in years.}
#'   \item{PLATINUM.STATUS}{The tumor's response to platinum therapy. Either
#'   "SENSITIVE" or "RESISTANT".}
#'   \item{TUMORRESIDUALDISEASE}{}
#'   \item{SURVIVALSTATUS}{The survival status of the patient. Either "LIVING"
#'   or "DECEASED".}
#'   \item{SURVIVALMONTHS}{}
#'   \item{RECURRENCE.MONTHS}{}
#'   \item{Differentiated.binary}{\code{TRUE} if \code{SUBTYPE} is
#'   "Differentiated"; \code{FALSE} otherwise.}
#'   \item{Immunoreactive.binary}{\code{TRUE} if \code{SUBTYPE} is
#'   "Immunoreactive"; \code{FALSE} otherwise.}
#'   \item{Mesenchymal.binary}{\code{TRUE} if \code{SUBTYPE} is
#'   "Mesenchymal"; \code{FALSE} otherwise.}
#'   \item{Proliferative.binary}{\code{TRUE} if \code{SUBTYPE} is
#'   "Proliferative"; \code{FALSE} otherwise.}
#'   \item{SILHOUETTE.WIDTH}{}
#'   \item{fake}{}
#' }
#'
#' @section fData:
#' \describe{
#'   \item{RefSeq}{The NCBI reference sequence.}
#' }
#'
#' @source CPTAC study
#'
#' @examples
#' data("cptac_oca")
#' head(pData(oca.set))
#' head(fData(oca.set))
#' head(exprs(oca.set))
"oca.set"


#' Example of a Longitudinal Biomarker Study
#'
#' Targeted selected reaction monitoring analysis. Time points A, B, C, and D
#' increase in the corresponding order. A & B are prior to the disease
#' diagnosis for the cases. C & D are correspondingly after the diagnosis.
#' Time point D is after the treatment.
#'
#' @docType data
#'
#' @usage data("longitudinal_biomarker_study")
#'
#' @format An \code{MSnSet} with 300 features and 236 samples.
#'
#' @section pData:
#' A data frame consisting of 236 rows (samples) and 12 columns:
#' \describe{
#'   \item{Sample}{The sample identifier. This is the \code{SubjID} and
#'   \code{TimePoint} concatenated with an underscore.}
#'   \item{isQC}{\code{TRUE} if the sample is a quality control sample.}
#'   \item{MatchID}{The unique subject identifier.}
#'   \item{Type}{"Control" or "Case".}
#'   \item{Plate}{The microplate identifier. Either "01", "02", "03", or "04".}
#'   \item{Well}{The microplate well location. \code{PlateRow} concatenated
#'   with \code{PlateCol}.}
#'   \item{PlateRow}{The row of the microplate where the sample is located.
#'   Letters A-H.}
#'   \item{PlateCol}{The column of the microplate where the sample is located.
#'   Strings "02" to "12".}
#'   \item{TimePoint}{Time point: A, B, C, or D.}
#'   \item{Days}{Number of days since disease diagnosis.}
#'   \item{SubjID}{Subject ID. This is the \code{MatchID} and a letter
#'   representing the \code{Type} ("C" for Control and "D" for Case)
#'   concatenated with an underscore.}
#'   \item{Age}{Subject age in years.}
#' }
#'
#' @section fData:
#' A data frame consisting of 300 rows (features) and 4 columns:
#' \describe{
#'   \item{Organism}{The organism name: "YEAST", "BOVIN", or "HUMAN".}
#'   \item{Protein}{The gene symbol and organism name concatenated
#'   with an underscore.}
#'   \item{Peptide}{The peptide sequence.}
#'   \item{isSpike}{\code{TRUE} if \code{Organism} is "YEAST" or "BOVIN".}
#' }
#'
#' @source xxx study
#'
#' @author Vlad Petyuk, 2018-11-09
#'
#' @examples
#' data("longitudinal_biomarker_study")
#' head(pData(longitudinal_biomarker_study))
#' head(fData(longitudinal_biomarker_study))
#' head(exprs(longitudinal_biomarker_study))
"longitudinal_biomarker_study"

