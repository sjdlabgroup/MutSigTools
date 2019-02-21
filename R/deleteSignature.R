#' @title Delete selected signatures from an existing set of mutation signatures.
#' 
#' @description Delete selected signatures from an existing set of mutation signatures. 
#' @usage deleteSignature(sigmatrix, del_sig)  
#' @param sigmatrix An object of class \emph{mutsig} describing the existing set of signatures.  
#' @param del_sig signature(s) to be removed.
#' @keywords delete mutation signatures. 
#' @return An object of \emph{mutsig} class with the deleted signatures.
#' @details Deletion of one or more mutation signatures from an existing set of mutation signatures upon which the updated signature set is returned.
#' @examples 
#' reducedSigMatrix=deleteSignature(sigmatrix=signatures.cosmic, c(2,5)) # delete signature 2 and 5.
#' @seealso \code{\link{addSignature}}, \code{\link{renameSignature}} and \code{\link{mergeSignature}}.
#' 
deleteSignature <-function(sigmatrix,del_sig) {
  return(t(subset(t(sigmatrix),select=-del_sig)))
}
