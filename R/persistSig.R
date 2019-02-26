#' @title Mutation signatures in different allele frequency spectrum
#' 
#' @description Determine the burden of different mutation signatures acrross different allele frequency ranges.
#'
#' @param snv A dataframe having \emph{sample, chr, pos, ref, alt and/or freq} as it's columns. This dataframe may have more than one sample.
#' @param th_vec_lw A vector of lower limits of frequency ranges.
#' @param th_vec_up A vector of upper limits of frequency ranges.
#' @return A list of samples having mutational signature corresponding to each allele frequency range
#' 
#' @export
#' @examples 
#' data(snv_sample)   # load 'snv' dataframe object
#' mut_sig_per_freq_range=persistSig(snv=snv_sample,th_vec_lw=c(0,0.4), th_vec_up=c(0.1,1))  # list of samples having mutational signature for each allele frequency range i.e., in our example (0.0 - 0.1) & (0.4 - 1.0).
#' 
#' @seealso \code{\link{vcfToSNV}} to generate \emph{snv} dataframe.

persistSig  <- function(snv,th_vec_lw,th_vec_up, bsg=NULL)
{
  ###### function for finding signature #######
  
  sig_find <- function(sample_name,sigs.input ){
    col_order=colnames(sigs.input)
    sample_sigs.input=sigs.input[sample_name,]
    decons_out=whichSignatures(tumor.ref = sample_sigs.input  , sample.id=sample_name, signatures.ref = signatures.cosmic, associated = c(), signatures.limit = NA, signature.cutoff = 0.00, contexts.needed = TRUE, tri.counts.method = "default")
    t0=as.vector(unlist(decons_out$weights, use.names=FALSE))
    names(t0)=c(paste0('sig.',seq(1:30))) #mut_burden_col_names
    return(t0)
  }
  ###### function to calculate mutation weight for all ranges of a sample #########
  all_freq_range_sig_list <- function(snv,sample_name,th_vec_lw,th_vec_up, bsg) {
    sample_sig_mat_at_ths=matrix(data=NA, nrow=length(th_vec_lw),ncol=30)
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
      sig_vec=sig_find(sample_name,sigs.input)
      sample_sig_mat_at_ths[th_range,]=sig_vec
      rm(sample_snv_freq_range)
    }
    colnames(sample_sig_mat_at_ths)=paste0('sig.',seq(1:30))
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
  mut_sig_list=lapply(samples, function(x) all_freq_range_sig_list(snv,x,th_vec_lw,th_vec_up,bsg))    # for each sample
  names(mut_sig_list)=samples
  return(mut_sig_list)
}
