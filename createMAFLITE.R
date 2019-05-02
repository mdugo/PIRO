

createMAFLITE <- function(fileName) {
  df <- read.table(
    fileName,
    header = T,
    sep = "\t",
    as.is = T,
    quote = ""
  )
  
  maflite <- data.frame(
    chr = df$CHROM,
    start = df$OPOS,
    end = df$OPOS,
    ref_allele = df$OREF,
    alt_allele = df$OALT,
    type = df$TYPE,
    stringsAsFactors = F
  )
  maflite$start[maflite$type == "ins"] <-maflite$start[maflite$type == "ins"] - 1
  maflite$end[maflite$type == "del"] <-maflite$start[maflite$type == "del"] + (nchar(maflite$ref_allele[maflite$type == "del"])-1)
  maflite$end[maflite$type == "mnp"]<-maflite$end[maflite$type == "mnp"] + (nchar(maflite$ref_allele[maflite$type == "mnp"])-1)
  maflite <- maflite[, colnames(maflite) != "type"]
  write.table(
    maflite,
    file = gsub("_vcf_tsv_merged_properview.txt", "_maflite.txt", fileName),
    sep = "\t",
    row.names = F,
    quote = F
  ) 
}