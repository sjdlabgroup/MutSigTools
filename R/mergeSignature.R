#' @title Creates a new mutational signature by taking a linear combination of other set of signatures.
#'
#' @description This function creates a new mutational signature by taking a linear combination of k other signatures s11,s12,...,skn located in a matrix of mutational signatures. In order to use this function, a weight vector w1, w2...wk must be specified. The new mutational signature is created by multiplying the component signature values by the  corresponding weights. The formula for the entries of the new signature is s11w1+s21w2+...+sk1wk, s12w1+s22w2+...+sk2wk,...s1nw1+s2nw2+...+sknwk. The function returns this new signature plus all of the signatures in the original dataset that were not merged into the new matrix of signatures.    
#' @usage mergeSignature(sigmatrix,sigs,weights)      
#' @param sigmatrix A \emph{mutsig} object of mutational signatures.
#' @param sig A vector that describes the signatures to be merged.
#' @param weights A vector that describes the relative weights of each of the component signatures in the merged mutational signature. 
#' @keywords merge mutation signatures.
#' @return A new matrix of signatures consisting of the signatures in the original matrix of mutational signatures that were not merged as well as the mutational signature that is a result of merging the input signatures. For example, if there were 10 signatures in the original signature matrix and signatures 1-6 were set to be merged, then the output value would be a matrix of 5 signatures(signature 7,8,9,10 and the result of merging signatures 1-6). 
#' @examples 
#' weights=rep(1/27,27)      
#' mergeSignature(signatures.cosmic,1:27,weights)      
#' weights=rep(1/6,6)  
#' mergeSignatures(signatures.cosmic,1:6,weights)   
#' 
#' @seealso \code{\link{addSignature}}, \code{\link{renameSignature}} and \code{\link{deleteSignature}}.

mergeSignature <-function(sigmatrix,sig,weights) {
  mergesig=colSums(sigmatrix[sig,]*weights)
  return(rbind(t(subset(t(sigmatrix),select=-sig)),mergesig))
}