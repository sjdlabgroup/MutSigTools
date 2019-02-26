#' @title Rename selected signatures from an existing set of signatures.    
#' 
#' @description Rename selected signatures from an existing signature-set.
#' @usage renameSignature(sigmatrix,selectSig,renameSig)
#' @param sigmatrix An object of class \emph{mutsig} describing the existing set of signatures.  
#' @param selectSig A vector containing name(s) of the signature(s) to be removed.
#' @param renameSig A vector containing updated name(s) of the signatures.
#' @keywords rename mutation signatures. 
#' @return An object of \emph{mutsig} class.
#' @export
#' @details Renaming of selected mutation signatures to an existing signature set, upon which the updated signature set is returned.
#' @examples 
#' renamedSigMatrix=renameSignature(sigmatrix=signatures.cosmic, selectSig=c(2,5), 
#' renameSig=c('two', 'five')) # rename signature 2 and 5.
#' @seealso \code{\link{addSignature}}, \code{\link{deleteSignature}} and \code{\link{mergeSignature}}.
#' 
renameSignature <-function(sigmatrix,selectSig,renameSig) {
  rownames(sigmatrix)[selectSig]=renameSig
  return(sigmatrix)
}
