#' @title Identify enrichment for signatures in individual samples relative to a set of reference samples.
#'
#' @description Determine over-represented mutation signatures in individual case sample(s) relative to that in the panel of control samples. This is suitable when case samples are heterogeneous, and only some of them might have excess of certain mutation signatures of interest relative to the control population..
#'
#' @usage enrichSig(contextfreq.cases, contextfreq.controls, signatures.ref=signatures.cosmic, threshold=0.05)
#' 
#' @param contextfreq.cases A data frame of class \emph{contextfreq} containing mutation frequency in tri-nucleotide contexts in case samples
#' @param contextfreq.controls A data frame of class \emph{contextfreq} containing mutation frequency in tri-nucleotide contexts in control samples
#' @param signatures.ref An object of class mutsig comprising the set of signatures. \emph{(signatures.nature2013 or signatures.cosmic or signatures.cosmic.2019 )}, Default: \emph{'signatures.cosmic'} 
#' @param threshold Threshold for uncorrected percentile score. Default: \emph{0.05} 

 
#' @keywords enriched signatures  
#' @return An object of class \emph{enrichSig.obj} providing the following information:
#' \enumerate{ 
#'   \item \strong{n.case}:            Number of case samples. 
#'   \item \strong{n.control}:		      Number of control samples. 
#'   \item \strong{case.weights}:		  A data frame containing estimated weights of known mutation signatures in the case samples. 
#'   \item \strong{control.weights}:		A data frame containing estimated weights of known mutation signatures in the control samples.
#'   \item \strong{case.percentile}:		Percentile scores corresponding to the extent of enrichment of known mutation signatures in case sample(s) relative to that in the control samples.
#'   }
#'   
#' @details Determine over-represented mutation signatures in individual case sample(s), highlighting those that are significantly enriched. The extent of enrichment is indicated using a percentile score, with low scores indicating high enrichment for specific mutation signatures in a case sample relative to that in the panel of control samples.
#' @export
#' @examples
#' data(contextfreq.cases_test)  # load case data
#' data(contextfreq.controls_test) # load control data
#' enrich_obj=enrichSig(contextfreq.case=contextfreq.cases_test, contextfreq.controls=contextfreq.controls_test, 
#' signatures.ref=signatures.cosmic, threshold=0.05)
#' @seealso \itemize{
#' \item \code{\link[deconstructSigs]{signatures.cosmic}}, \code{\link{confidenceSig}} for robust signatures and \code{\link{caseControlSig}} to identify signatures with significantly higher mutation burden in case samples over control samples. 
#' \item To generate \emph{contextfreq} object from \emph{snv} dataframe use \code{\link{processSNV}} and \code{\link{vcfToSNV}}.
#' }


enrichSig <-function(contextfreq.cases, contextfreq.controls, signatures.ref=signatures.cosmic, threshold=0.05) {
  control.weights <-c()
  for(i in 1:nrow(contextfreq.controls)) {
    control.weights=rbind(control.weights,whichSignatures(contextfreq.controls[i,],sample.id = rownames(contextfreq.controls[i,]), signatures.ref, contexts.needed = T)$weights)
  }
  case.weights <-c()
  case.p <-c()
  for(i in 1:nrow(contextfreq.cases)) {
    case.weights <-rbind(case.weights,whichSignatures(contextfreq.cases[i,],sample.id = rownames(contextfreq.cases[i,]), signatures.ref, contexts.needed = T)$weights)
    case.p <-rbind(case.p, 1-apply(rbind(case.weights[i,],control.weights),2,rank)[1,]/(nrow(contextfreq.controls)+1))
  }
  case.p[case.p > threshold] <- "."
  
  enrichSig.obj <-list()
  enrichSig.obj$n.case <-nrow(contextfreq.cases)
  enrichSig.obj$n.control <-nrow(contextfreq.controls)
  enrichSig.obj$case.weights <-data.frame(case.weights)
  enrichSig.obj$control.weights <-data.frame(control.weights)
  enrichSig.obj$case.percentile <-data.frame(case.p)
  rownames(enrichSig.obj$case.percentile) <-rownames(contextfreq.cases)
  return(enrichSig.obj)
}
