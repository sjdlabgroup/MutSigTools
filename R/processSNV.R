#' @title Estimate mutation frequency in tri-nucleotide motif contexts.
#'
#' @description Estimate frequency of point mutations in tri-nucleotide motif contexts for one or more input samples.
#' @usage processSNV(snv, bsg = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
#' @param snv A dataframe having \emph{sample, chr, pos, ref, alt and/or freq} as its columns. The \emph{snv} dataframe can be created by the \code{\link{vcfToSNV}} function.
#' @param bsg An object of class \emph{BSGenome}. Default: \emph{BSgenome.Hsapiens.UCSC.hg19::Hsapiens} 
#' 

#' @keywords mutation catalogue
#' @return A data frame of class \emph{contextfreq} containing mutation frequency in user-specified nucleotide contexts.
#' @details Estimation of mutation frequency in pre-determined, user-specified nucleotide motif contexts for input samples. The BSgenome object contains the reference genome information, and should be specified for alternative assembly of human genomes or nonhuman reference genomes. The output returns a data frame containing estimated mutation frequencies in the different genomic contexts, as provided in the \emph{contextfreq} object, for each input sample.
#' 
#' @export
#' 
#' @examples
#' load(file = "data/snv_sample.rda")
#' context.freq=processSNV(snv=snv_sample, bsg = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' @seealso  \code{\link{vcfToSNV}} to generate \emph{snv} dataframe, \code{\link[BSgenome]{BSgenome}} and  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}. 

processSNV <- function(snv, bsg = BSgenome.Hsapiens.UCSC.hg19::Hsapiens){
  
  if(ncol(snv)<5)
  {
    t=c("sample","chr","pos","ref","alt")
    col_absent=setdiff(t,colnames(snv))
    err_message=paste0("INPUT snv is missing '",paste(col_absent,collapse = "', '"),"' cloumns")
    stop(err_message, call. = TRUE, domain = NULL)
    geterrmessage()
    
  }
  
snv <- snv[which(snv$ref %in% c('A', 'T', 'C', 'G') & snv$alt %in% c('A', 'T', 'C', 'G')),]

snv$chr <- factor(snv$chr)
levels(snv$chr) <- sub("^([0-9XY])", "chr\\1", levels(snv$chr))
levels(snv$chr) <- sub("^MT", "chrM", levels(snv$chr))
levels(snv$chr) <- sub("^(GL[0-9]+).[0-9]", "chrUn_\\L\\1", levels(snv$chr), perl = T)

# Check the genome version the user wants to use
# If set to default, carry on happily
if(is.null(bsg)){
      genome.ref=BSgenome.Hsapiens.UCSC.hg19::Hsapiens
}
if(class(bsg) != 'BSgenome'){
      stop('The reference genome does not match with a BSgenome object.')
} else genome.ref=bsg

# Exclude chromosomes that do not exist in the BSGenome object
unknown.regions <- levels(snv$chr)[which(!(levels(snv$chr) %in% GenomeInfoDb::seqnames(genome.ref)))]
if (length(unknown.regions) > 0) {
      unknown.regions <- paste(unknown.regions, collapse = ',\ ')
      warning(paste('Data includes chromosomes that do not exist in the BSGenome object\n', unknown.regions, sep = ' '))
      snv <- snv[snv$chr %in% GenomeInfoDb::seqnames(genome.ref), ]
}
if(!all(snv$ref %in% c('A', 'T', 'C', 'G') & snv$alt %in% c('A', 'T', 'C', 'G'))) {
      stop("Only SNVs with standard base substitutions are supported.")
}
# Notify if the reference base does not match with that in the reference genome
refN = BSgenome::getSeq(genome.ref, snv$chr, snv$pos, snv$pos, as.character = T)
if(any(refN != snv$ref)) {
      warning(sprintf("Discordant reference base call in %d cases", sum(refN != snv$ref)))
}

# Extract trinucleotide contex##
contextfreq.sample=mut.to.sigs.input(mut.ref = snv, sample.id = "sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt")
return(contextfreq.sample)


}

