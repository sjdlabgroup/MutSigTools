#' @title Mutation signatures in different allele frequency spectrum
#' 
#' @description Determine the burden of different mutation signatures acrross different allele frequency ranges.The function should be used with caution; Local copy number and purity-corrected variant allele frequencies should be provided in the column 6 of the \emph{snv}
#'
#' @param snv A dataframe having \emph{sample, chr, pos, ref, alt and/or freq} as it's columns. This dataframe may have more than one sample.
#' @param th_vec_lw A vector of lower limits of frequency ranges.
#' @param th_vec_up A vector of upper limits of frequency ranges.
#' @param signatures.ref An object of class \emph{mutsig} comprising the set of signatures. \emph{(signatures.nature2013 or signatures.cosmic or signatures.cosmic.2019 )}, \emph{Default: 'signatures.cosmic'} 
#' @param bsg An object of class \emph{BSGenome}. Default: \emph{BSgenome.Hsapiens.UCSC.hg19} 

#' @return A list of samples having mutational signature corresponding to each allele frequency range
#' 
#' @export
#' @examples 
#' data(snv_sample)   # load 'snv' dataframe object
#' mut_sig_per_freq_range=persistSig(snv=snv_sample,th_vec_lw=c(0,0.4), th_vec_up=c(0.1,1), signatures.ref = signatures.cosmic,  bsg=BSgenome.Hsapiens.UCSC.hg19::Hsapiens)  # list of samples having mutational signature for each allele frequency range i.e., in our example (0.0 - 0.1) & (0.4 - 1.0).
#' 
#' @seealso \code{\link{vcfToSNV}} to generate \emph{snv} dataframe.

persistSig  <- function(snv,th_vec_lw,th_vec_up, signatures.ref = signatures.cosmic, bsg=BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
{
  ###### function for finding signature #######
  sigmat=signatures.ref
  sig_find <- function(sample_name,sigs.input,bsg, sigmat ){
    col_order=colnames(sigs.input)
    sample_sigs.input=sigs.input[sample_name,]
    decons_out=whichSignatures(tumor.ref = sample_sigs.input  , sample.id=sample_name, signatures.ref = sigmat, associated = c(), signatures.limit = NA, signature.cutoff = 0.00, contexts.needed = TRUE, tri.counts.method = "default")
    t0=as.vector(unlist(decons_out$weights, use.names=FALSE))
    names(t0)=rownames(sigmat)
    #names(t0)=c(paste0('sig.',sig_names)) #mut_burden_col_names
    return(t0)
  }
  ###### function to calculate mutation weight for all ranges of a sample #########
  all_freq_range_sig_list <- function(snv,sample_name,th_vec_lw,th_vec_up, bsg, sigmat) {
    sample_sig_mat_at_ths=matrix(data=NA, nrow=length(th_vec_lw),ncol=nrow(sigmat))
    rownames(sample_sig_mat_at_ths)=paste0('row',seq(1,length(th_vec_lw)))
    for (i in seq(1,length(th_vec_lw))){
      th_range=paste0(th_vec_lw[i],'-',th_vec_up[i])
      rownames(sample_sig_mat_at_ths)[i]=th_range
      sample_snv_freq_range=subset(snv,freq > th_vec_lw[i] & freq < th_vec_up[i] & sample== sample_name)
      if(isEmpty(sample_snv_freq_range$sample)==1){
        print(paste0(sample_name,': sample do nat have any mutation in ', th_range ,'allele frequency range'))
        next
      }
      sigs.input = mut.to.sigs.input(mut.ref = sample_snv_freq_range[,c(1,2,3,4,5)], sample.id = "sample",chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = bsg)
      sig_vec=sig_find(sample_name,sigs.input,bsg, sigmat)
      sample_sig_mat_at_ths[th_range,]=sig_vec
      rm(sample_snv_freq_range)
    }
    #colnames(sample_sig_mat_at_ths)=paste0('sig.',seq(1:30))
    colnames(sample_sig_mat_at_ths)=rownames(sigmat)
    
    return(sample_sig_mat_at_ths)
  }
  ###########  Main Function ########
  if(ncol(snv)!=6)
  {
    t=c("sample","chr","pos","ref","alt","freq")
    col_absent=setdiff(t,colnames(snv))
    err_message=paste0("INPUT snv is missing '",paste(col_absent,collapse = "', '"),"' cloumns")
    stop(err_message, call. = TRUE, domain = NULL)
    geterrmessage()

  }
  samples=as.vector(unique(snv$sample))
  mut_sig_list=lapply(samples, function(x) all_freq_range_sig_list(snv,x,th_vec_lw,th_vec_up,bsg, sigmat))    # for each sample
  names(mut_sig_list)=samples
  return(mut_sig_list)
}
