#' @title Create three PCA plots for signatures.
#'
#' @description Constructs three PCA plots showing variations among the signatures in terms of the weights of the principal components for different mutation signatures, as well as the eigenvalues and cosine similarity among expected and actual signatures.
#'
#' @usage signaturePCA(sigmat, pngfile )
#' @param sigmat An object of class \emph{mutsig} describing a set of signatures.
#' @param pngfile name of \emph{png} formatted image file to be created
#' @keywords mutation signamture PCA
#' @return Create three PCA plots with the following extensions to show difference among mutational signatures.
#' \enumerate{  
#'  \item \strong{*_Eigen.png} The eigenvalues of the PCA and the relative amount of the variation each eigenvalue.
#'  \item \strong{ *_var.png} The different single-nucleotide alteration components in terms of the two dimensions.
#'  \item \strong{*_PCA.png} The different signatures in terms of the 2 non-dimensional vectors derived by the PCA analysis and the cosine similarity between the estimate of each mutation signature using PCA and the actual mutation signature(in order to determine how accurately the 2 dimensions calculated by PCA represent the mutational signatures). 
#'  }   
#' @export
#' @examples
#' sigmat=signatures.cosmic
#' signaturePCA(sigmat, pngfile="signaturePCA")
#' 
#' @seealso \code{\link{signatureHeatmap}}, \code{\link[deconstructSigs]{signatures.cosmic}}, \code{\link{deleteSignature}}, \code{\link{renameSignature}} and \code{\link{mergeSignature}}.

signaturePCA <-function(sigmat, pngfile="signature") {
  sig.pca <- prcomp(sigmat,center = TRUE, scale. = TRUE)
  png(filename=paste0(pngfile,"_PCA.png"), height=600, width=600, res=300)
  fviz_pca_ind(sig.pca,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  dev.off()
  png(filename=paste0(pngfile,"_Var.png"), height=600, width=600, res=300)
  fviz_pca_var(sig.pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  dev.off()
  png(filename=paste0(pngfile,"_Eigen.png"), height=600, width=600)
  fviz_eig(sig.pca)
  dev.off()
  return(sig.pca)
}