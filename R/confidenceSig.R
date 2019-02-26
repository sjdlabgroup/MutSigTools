#' @title Provides interval of uncertainty for estimated weights of mutation signatures
#'
#' @description Provides an interval of uncertainty for estimated weights of known mutation signatures in the catalog of mutations from a set of samples.
#'
#' @usage confidenceSig(contextfreq.sample, subsample=0.8, iter=1000 , signatures.ref=signatures.cosmic, lbound=0.1, ubound=0.9, replace=FALSE)
#' 
#' @param contextfreq.sample A sample from the dataframe of class \emph{contextfreq} containing mutation frequency in trinucleotide contexts.
#' @param subsample Proportion of mutations included during each subsampling. Default: \emph{0.8}  (80 percent)
#' @param iter Number of iterations of subsampling. Default: \emph{1000}
#' @param replace should sampling be with replacement? \emph{TRUE or FALSE}
#' 
#' @param signatures.ref An object of class \emph{mutsig} comprising the set of signatures. Default: \emph{'signatures.cosmic'} 
#' @param lbound Lower bound of the interval of uncertainty for estimated weights of the signatures Default: \emph{0.1} (10 percent)
#' @param ubound Upper bound of the interval of uncertainty for estimated weights of the signatures Default: \emph{0.9} (90 percent)
#' @keywords uncertainty for estimated weights of mutation signatures 
#' @return An object containing the following information:
#' \enumerate{ 
#'   \item \strong{observed.weights:} A data frame containing estimated weights of known mutation signatures based on all mutations in a sample.
#'   \item \strong{median.weights:} A data frame containing estimated median of the weights of known mutation signatures based on iteratively subsampled mutations.
#'   \item \strong{ubound.weights:} A data frame containing upper-bound values of the weights of known mutation signatures based on iteratively subsampled mutations.
#'   \item \strong{lbound.weights:} A data frame containing lower-bound values of the weights of known mutation signatures based on iteratively subsampled mutations.
#'   }
#' @export
#' @details Provides an interval of uncertainty for estimated weights of known mutation signatures in the catalog of mutations from a set of samples. First, based on the catalog of mutations in a sample, weights of the mutation signatures are estimated. Next, mutations are subsampled iteratively without replacement, each time estimating the weights of the mutation signatures. The intervals of uncertainty for weights of each mutation signatures are determined by aggregating observations from a given number of iterations. 
#' @examples
#' data(contextfreq.sample_test)
#
#' robust_sig_object=confidenceSig(contextfreq.sample=contextfreq.sample_test[1,], subsample=0.8, 
#' iter=50, signatures.ref=signatures.cosmic, lbound=0.1, ubound=0.9, replace=FALSE)
#' @seealso 
#' \itemize{ 
#' \item \code{\link[deconstructSigs]{signatures.cosmic}}, \code{\link{caseControlSig}} to identify signatures with significantly higher mutation burden in case samples over control samples and \code{\link{enrichSig}} for enriched signatures in individual case samples.  
#'  \item To generate \emph{contextfreq} object from \emph{snv} dataframe use \code{\link{processSNV}} and \code{\link{vcfToSNV}}.
#'  }


confidenceSig <-function(contextfreq.sample, subsample=0.8, iter=1000, signatures.ref=signatures.cosmic, lbound=0.1, ubound=0.9, replace=FALSE) {
  mutevents <- rep(colnames(contextfreq.sample), contextfreq.sample[1,])
  permwt=c()
  for(i in 1:iter) {
    mutsample <- mutevents[sample(1:length(mutevents), round(length(mutevents)*subsample), replace=replace)]
    table(factor(mutsample, levels=colnames(contextfreq.sample)))
    
    contextfreq.subsample <-data.frame(t(as.numeric(table(factor(mutsample, levels=colnames(contextfreq.sample))))))
    colnames(contextfreq.subsample) <-colnames(contextfreq.sample)
    rownames(contextfreq.subsample) <- i
    permwt=rbind(permwt,whichSignatures(contextfreq.subsample,sample.id = i, signatures.ref, contexts.needed = T)$weights)
  }
  mutsig.obj=list()
  mutsig.obj$observed.weights <- whichSignatures(contextfreq.subsample, sample.id = i, signatures.ref, contexts.needed = T)$weights
  mutsig.obj$subsampling <- paste0(subsample*100,"%")
  mutsig.obj$iter <- iter
  mutsig.obj$confidence.interval<-paste0("(",ubound,",",lbound,")")
  mutsig.obj$median.weights <- apply(contextfreq.sample,2,function(x){ quantile(x, probs=c(0,0.5,1))[2]})
  mutsig.obj$ubound.weights <-apply(contextfreq.sample,2,function(x){ quantile(x, probs=c(0,ubound,1))[2]})
  mutsig.obj$lbound.weights <-apply(contextfreq.sample,2,function(x){ quantile(x, probs=c(0,lbound,1))[2]})
  return(mutsig.obj)
}
