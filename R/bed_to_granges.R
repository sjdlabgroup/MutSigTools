
bed_to_granges <- function(file){
  df <- read.table(file, header=F,stringsAsFactors=F)

  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }

  if(length(df)<3){
    stop("File has less than 3 columns")
  }

  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]

  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }


  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}
