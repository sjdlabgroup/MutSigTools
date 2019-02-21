#' @title Identify signatures with significantly higher mutation burden in cases over controls.
#' 
#' @description Identify signatures with significantly higher mutation burden in case samples over control samples.
#'
#' @usage caseControlSig(contextfreq.cases, contextfreq.controls, signatures.ref=signatures.cosmic, threshold=0.05, adjust="fdr")
#' 
#' @param contextfreq.cases A data frame of class \emph{contextfreq} containing mutation frequency in tri-nucleotide contexts in case samples
#' @param contextfreq.controls A data frame of class \emph{contextfreq} containing mutation frequency in tri-nucleotide contexts in control samples
#' @param signatures.ref An object of class \emph{mutsig} comprising the set of signatures. \emph{Default: 'signatures.cosmic'} 
#' @param threshold Threshold for uncorrected percentile score. \emph{Default: 0.05}
#' @param adjust Method for p-value correction for multiple testing. Options are as provided in the function \emph{p.adjust}. The default method is the \emph{"FDR"} method.  

#' @keywords enriched signatures  
#' @return An object providing the following information:
#' \enumerate{ 
#'   \item \strong{n.case}:             Number of case samples. 
#'   \item \strong{n.control}:		      Number of control samples. 
#'   \item \strong{case.weights}:		    A data frame containing estimated weights of known mutation signatures in the case samples. 
#'   \item \strong{control.weights}:		A data frame containing estimated weights of known mutation signatures in the control samples.
#'   \item \strong{p.value}:		        Indicates statistical significance of higher mutation burden of signatures in cases over controls.
#'   \item \strong{adjust.p.value:}		  P-values adjusted for multiple testing correction.
#'   }
#'   
#' @details Identifes signatures with significantly higher mutation burden in case samples campared to control samples using the Wilcoxon Rank sum test. This is suitable when the samples are homogeneous within respective groups, and signatures significantly over-represented in cases relative to the controls are of interest.
#' @export
#' @examples
#' load(file = "data/contextfreq.cases_test.rda")
#' load(file = "data/contextfreq.controls_test.rda")
#' signif_signatures_obj=caseControlSig(contextfreq.cases=contextfreq.cases_test, contextfreq.controls=contextfreq.controls_test, signatures.ref=signatures.cosmic, threshold=0.05, adjust="fdr")
#' @seealso 
#' \itemize{ 
#' \item \code{\link[deconstructSigs]{signatures.cosmic}}, \code{\link{confidenceSig}} for robust signatures and \code{\link{enrichSig}} for enriched signatures in individual case samples.
#' \item To generate \emph{contextfreq} object from  \emph{snv} dataframe use \code{\link{processSNV}} and \code{\link{vcfToSNV}}.
#' }

caseControlSig <-function(contextfreq.cases, contextfreq.controls, signatures.ref=signatures.cosmic, threshold=0.05, adjust="fdr") {
  control.weights <-c()
  for(i in 1:nrow(contextfreq.controls)) {
    control.weights=rbind(control.weights,whichSignatures(contextfreq.controls[i,],sample.id = rownames(contextfreq.controls[i,]), signatures.ref, contexts.needed = T)$weights)
  }
  
  case.weights <-c()
  for(i in 1:nrow(contextfreq.cases)) {
    case.weights <-rbind(case.weights,whichSignatures(contextfreq.cases[i,],sample.id = rownames(contextfreq.cases[i,]), signatures.ref, contexts.needed = T)$weights)
  }
  
  caseControl.p <-c()
  caseControl.adjust.p <-c()
  for(i in 1:nrow(signatures.ref)) {
    caseControl.p <- c(caseControl.p,wilcox.test(t(case.weights[i]),t(control.weights[i]), alternative = "g")$p.value)
  }
  caseControl.adjust.p <-p.adjust(caseControl.p, method=adjust)
  caseControl.p <- as.numeric(formatC(caseControl.p, format="e", digits=3))
  caseControl.adjust.p <- as.numeric(formatC(caseControl.adjust.p, format="e", digits=3))
  caseControl.p[caseControl.p > threshold] <- "."
  caseControl.adjust.p[caseControl.adjust.p > threshold] <- "."
  
  caseControlSig.obj <-list()
  caseControlSig.obj$n.case <-nrow(contextfreq.cases)
  caseControlSig.obj$n.control <-nrow(contextfreq.controls)
  caseControlSig.obj$case.weights <-data.frame(case.weights)
  caseControlSig.obj$control.weights <-data.frame(control.weights)
  caseControlSig.obj$p.value <-data.frame(t(caseControl.p))
  colnames(caseControlSig.obj$p.value) <- rownames(signatures.ref)
  caseControlSig.obj$adjust.p.value <-data.frame(t(caseControl.adjust.p))
  colnames(caseControlSig.obj$adjust.p.value) <- rownames(signatures.ref)
  caseControlSig.obj$adjust.method <-adjust
  return(caseControlSig.obj)
}
