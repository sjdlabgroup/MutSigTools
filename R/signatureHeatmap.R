#' @title Create a heatmap of motif weights for signatures.
#'
#' @description Construct a heatmap showing relative weights of tri-nucleotide motif for different mutation signatures.
#'
#' @usage signatureHeatmap(sigmat, pngfile=..., )
#' @param sigmat An object of class \emph{mutsig} describing the existing set of signatures.
#' @param pngfile The name of \emph{png} formatted image file to be created.
#' @keywords mutation signamture heatmap
#' @return A \emph{png} formatted image file.  
#' @export
#' @Details Create a heatmap to show relative weights of trinucleotide motifs for different mutation signatures in a signature-set. It is possible to use alternative motifs such as penta-nucleotide motifs etc.
#' @examples
#' sigmat=signatures.cosmic
#' signatureHeatmap(sigmat, pngfile="signature.png", mar=c(6,8)
#' @seealso \code{\link{signaturePCA}}, \code{\link{clusterSignature}} , \code{\link{addSignature}} and \code{\link[deconstructSigs]{signatures.cosmic}}


signatureHeatmap <-function(sigmat, pngfile="signature", mar=c(6,8)) {
  my_palette <- colorRampPalette(c("white", "pink", "red"))(n = 299)
  col_breaks = c(seq(0,0.02,length=100),
                 seq(0.03,0.3,length=100),
                 seq(0.31,1,length=100))
  png(filename=paste0(pngfile,'_heatmap.png'), height=600, width=1200)
  heatmap.2(as.matrix(sigmat), trace="none", key=F, scale="none", col=my_palette, breaks=col_breaks, margins=mar)
  dev.off()
}