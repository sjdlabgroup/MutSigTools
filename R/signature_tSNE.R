
#' @title Create tSNE plot of signatures.
#'
#' @description Constructs tSNE plots of signatures. 
#'
#' @usage signature_tSNE(sigmat, pngfile )
#' @param sigmat An object of class \emph{mutsig} describing a set of signatures.
#' @param pngfile name of \emph{png} formatted image file to be created
#' @keywords mutation signature tSNE
#' @return Create tSNE plot to show similarity among mutational signatures.
#' \enumerate{  
#'  \item \strong{*_tSNE.png} The tSNE plot of mutation signatures.
#'  }   
#' @export
#' @examples
#' sigmat=signatures.cosmic
#' signature_tSNE(sigmat, pngfile="signature_tSNE")
#' 
#' @seealso \code{\link{signatureHeatmap}}, \code{\link{signaturePCA}}, \code{\link[deconstructSigs]{signatures.cosmic}}, \code{\link{deleteSignature}}, \code{\link{renameSignature}} and \code{\link{mergeSignature}}.

signature_tSNE <-function(sigmat, pngfile="signature") {
  set.seed(20)
  sig_names=sapply(rownames(sigmat), function(x) strsplit(x,"[.]")[[1]][[2]])
  perpex=round((length(rownames(sigmat))/3)-1)
  
  png(filename=paste0(pngfile,"_tSNE.png"))
  tsne <- Rtsne(as.matrix(sigmat), dims = 2, perplexity=perpex, verbose=TRUE, max_iter = 500)
  p <- plot(tsne$Y, pch=20,  lwd = .01, col=gray(1), xlab='tSNE1',ylab='tSNE2')
  text(tsne$Y, labels=sig_names)
  print(p)

  dev.off()
  
}