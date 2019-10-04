#' @title Convert vcf file to snv dataframe
#'
#' @description Read vcf formatted file to catalog mutations in a data frame of class snv for downstream analysis. A snv dataframe consists of sample, chr, pos, ref, alt and/or freq columns.
#' 
#' @usage vcfToSNV(vcf, allelefreq = FALSE)  
#' @param vcf  The path of a standard vcf file.
#' @param allelefreq A logical input, TRUE if mutation or variant allele frequency is required in output object. Default : FALSE
#' @return A snv dataframe having datframe with column names as sample, chr, pos, ref, alt, freq. Where, filename will be considered as sample name. \strong{Note}: vcf file name with no extension is the sample name.

#' @details Read vcf formatted file to catalog mutations in a data frame of class snv for downstream analysis foreach input sample.
#' @export
#' @examples 
#' vcf_file=system.file("extdata", "test.vcf", package = "MutSigTools", mustWork = TRUE)
#' 
#' snv=vcfToSNV(vcf=vcf_file, allelefreq=TRUE)
#' 
#' head(snv)
#' 
#' @seealso \href{https://faculty.washington.edu/browning/intro-to-vcf.html}{VCF} file format.
#' 
#' 
vcfToSNV <- function(vcf, allelefreq=FALSE)
{
      
      vcf_data=read.vcf(vcf)

            info_split=strsplit(x=vcf_data$vcf$INFO,split = ';')
            sample_name_vec=rep(sub('\\.vcf$', '', basename(vcf)),nrow(vcf_data$vcf) )
            AF_exists_list=sapply(info_split, function(x) x[grep('AF',x)])

            if(allelefreq==FALSE) {
            	snv=data.frame(sample=sample_name_vec,chr=vcf_data$vcf$CHROM,pos=vcf_data$vcf$POS, ref=vcf_data$vcf$REF, alt=vcf_data$vcf$ALT)
            }
            else if(length(AF_exists_list[lapply(AF_exists_list,length)>0]) < nrow(vcf_data$vcf))
            {
                  warning("Allele frequency is missing for some mutations or formatted differently" )
                  snv=data.frame(sample=sample_name_vec,chr=vcf_data$vcf$CHROM,pos=vcf_data$vcf$POS, ref=vcf_data$vcf$REF, alt=vcf_data$vcf$ALT)

            } else {
              freq_vec=sapply(info_split, function(x) as.numeric(strsplit(x[grep('AF',x)],'=')[[1]][2]))
            	snv=data.frame(sample=sample_name_vec,chr=vcf_data$vcf$CHROM,pos=vcf_data$vcf$POS, ref=vcf_data$vcf$REF, alt=vcf_data$vcf$ALT, freq=freq_vec)
            }

      return(snv)
}
