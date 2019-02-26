#' @title Add new signature(s) to an existing set of mutation signatures
#' 
#' @description Add new signature(s) to an existing set of mutation signatures, and  return the updated signature set.
#' @usage addSignature(sigmatrix2, sigmatrix1)  
#' @param sigmatrix2 An object of class \emph{mutsig} describing the new signature.   
#' @param sigmatrix1 An object of class \emph{mutsig} describing the existing set of signatures.
#' @keywords Add mutation signatures to a set of signatures. 
#' @return An object of \emph{mutsig} class.
#' @export
#' @details Addition of one or more new mutation signature(s) to an existing set of mutation signatures upon which the updated signature set is returned.
#' @examples 
#' newSigMatrix=addSignature(sigmatrix2=signatures.cosmic[c('Signature.10','Signature.11'),], 
#' sigmatrix1=signatures.cosmic[c('Signature.1','Signature.2'),])
#' @seealso \code{\link{deleteSignature}}, \code{\link{renameSignature}} and \code{\link{mergeSignature}}.

addSignature <-function(sigmatrix2,sigmatrix1) {
  return(rbind(sigmatrix1, sigmatrix2))
}