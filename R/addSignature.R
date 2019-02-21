#' @title Add new signature(s) to an existing set of mutation signatures
#' 
#' @description Add new signature(s) to an existing set of mutation signatures, and  return the updated signature set.
#' @usage addSignature(sigmatrix1, sigmatrix2)  
#' @param sigmatrix2 An object of class \emph{mutsig} describing the new signature.   
#' @param sigmatrix1 An object of class \emph{mutsig} describing the existing set of signatures.
#' @keywords Add mutation signatures to a set of signatures. 
#' @return An object of \emph{mutsig} class.
#' @export
#' @details Addition of one or more new mutation signature(s) to an existing set of mutation signatures upon which the updated signature set is returned.
#' @examples 
#' addSignature(sigmatrix2, sigmatrix1)
#' @seealso \code{\link{deleteSignature}}, \code{\link{renameSignature}} and \code{\link{mergeSignature}}.

addSignature <-function(sigmatrix2,sigmatrix1) {
  return(rbind(sigmatrix2,sigmatrix1))
}