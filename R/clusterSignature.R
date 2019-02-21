#' Hierarchical clustering of mutational signature.
#' 
#' Performs a hierarchical clustering analysis of an input mutational signature matrix using the Euclidean distance metric.    
#' @usage clusterSignatures(sigmatrix)  
#' @param sigmatrix An input matrix of mutational signatures.  
#' @keywords cluster mutation signatures. 
#' @return A hierarchical clustering object formed from the input matrix of mutational signatures. This can then be plotted to show how the mutational signatures in the input matrix are related to one another.   
#' @examples 
#' hclust_obj=clusterSignature(sigmatrix=signatures.cosmic) 
#' \code{\link{signaturePCA}} and \code{\link{signaturesHeatmap}}

clusterSignature <-function(sigmatrix) {
  hc=hclust(dist(sigmatrix))
  return(hc)
}
