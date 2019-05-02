

###################################################################
# This function merges the vcf and tsv files obtained from a
# tumor-normal workflow of Ion Reporter.
# The function splits the information contained in the INFO and 
# SAMPLES columns of the vcf and merges the vcf with the tsv.
# The merged output is written in a tab-delimited text file.
# Redundant columns contained both in the vcf and tsv files are
# removed and only those contained in the vcf are kept.
# The effective position, reference and alternate allele for the 
# variant are in the OPOS, OREF and OALT columns.
#
# WARNING!!!!! This function works only for vcf and tsv files that 
# DO NOT CONTAIN CNV calls
###################################################################

############################ USAGE ################################

# mergeVcfTsv(sampleName,path.in,path.out,workflow)

# Arguments

# sampleName: character string for the sample name. The sample name
# must be contained in the vcf and tsv file names. Ideally is the 
# sample barcode

# path.in: the path to the folder containing the vcf and tsv files

# path.out: the path where the output will be written

######################### FUNCTION ################################

mergeVcfTsv<-function(sampleName,path.in,path.out,workflow=c("T-N","SS")){
  if(workflow=="T-N"){
  vcfFile<-grep(sampleName,list.files(path.in,pattern="Filtered",full.names=T,recursive=T),value=T)
  tsvFile<-grep(sampleName,list.files(path.in,pattern="\\.tsv",full.names=T,recursive=T),value=T)
  
  vcfTemp<-read.table(vcfFile,sep="\t",as.is=T,comment.char="#",header=F,stringsAsFactors=F)
  colnames(vcfTemp)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","TUMOR","NORMAL")
  vcfTemp<-vcfTemp[vcfTemp$ALT!="<CNV>",]
  
  tsvTemp<-read.table(tsvFile,sep="\t",header=T,comment.char="",skip=2,stringsAsFactors = F,quote="")
  tsvTemp<-tsvTemp[tsvTemp$type!="CNV",]
  
  if(nrow(vcfTemp)!=nrow(tsvTemp)){
    stop("Error: different number of variants in vcf and tsv files")
  }
  
  # splitting of INFO field -------------------------------------------------
  
  uniqueField<-unique(unlist(lapply(sapply(vcfTemp$INFO,strsplit,";"),function(x)gsub("=.*","",x))))
  infoTemp<-strsplit(vcfTemp$INFO,";")
  list2df<-function(x){
    y<-data.frame(field=gsub("=.*","",x),value=gsub(".*=","",x),stringsAsFactors = F)
  }
  infoTemp<-lapply(infoTemp,list2df)
  
  for(k in 1:length(infoTemp)){
    field=uniqueField[!uniqueField%in%infoTemp[[k]][,1]]
    if(length(field)>0){
      missing<-data.frame(field=field,value="")
      infoTemp[[k]]<-rbind(infoTemp[[k]],missing)
    }
    infoTemp[[k]]<-infoTemp[[k]][order(infoTemp[[k]][,1]),]
    infoTemp[[k]]<-infoTemp[[k]][infoTemp[[k]][,1]!="FUNC",]
  }
  
  if(length(table(sapply(infoTemp,nrow)))>1){
    stop("Error: the number of field in the INFO column is different across variants")
  } 
  columnNames<-as.vector(infoTemp[[1]]$field)
  infoTemp<-lapply(infoTemp,function(x)as.vector(x$value))
  dfInfo <- data.frame(matrix(unlist(infoTemp), nrow=length(infoTemp), byrow=T),stringsAsFactors = F)
  colnames(dfInfo)<-columnNames
  
  # attach to vcf data frame
  if(nrow(vcfTemp)!=nrow(dfInfo)){
    stop("Error: different number of rows between vcf and INFO data frame")
  }
  vcfTemp<-cbind(vcfTemp,dfInfo)

  # splitting of TUMOR field -------------------------------------------------
  
  tumorTemp<-lapply(strsplit(vcfTemp$TUMOR,":"),"[",1:2)
  dfTumor<-data.frame(GT_TUMOR=sapply(tumorTemp,"[",1),GQ_TUMOR=sapply(tumorTemp,"[",2),stringsAsFactors = F)

  # splitting of NORMAL field -------------------------------------------------
  
  normalTemp<-lapply(strsplit(vcfTemp$NORMAL,":"),"[",c(1,2,4,8,9))
  dfNormal<-data.frame(GT_NORMAL=sapply(normalTemp,"[",1),
                       GQ_NORMAL=sapply(normalTemp,"[",2),
                       FDP_NORMAL=sapply(normalTemp,"[",3),
                       FAO_NORMAL=sapply(normalTemp,"[",4),
                       AF_NORMAL=sapply(normalTemp,"[",5),
                       stringsAsFactors = F)
  
  # attach to vcf data frame
  if(nrow(vcfTemp)!=nrow(dfTumor)){
    stop("Error: different number of rows between vcf and TUMOR data frame")
  }
  vcfTemp<-cbind(vcfTemp,dfTumor)
  
  if(nrow(vcfTemp)!=nrow(dfNormal)){
    stop("Error: different number of rows between vcf and NORMAl data frame")
  }
  vcfTemp<-cbind(vcfTemp,dfNormal)
  
  vcfTemp<-vcfTemp[,!colnames(vcfTemp)%in%c("INFO","FORMAT","TUMOR","NORMAL")]
  
  # remove multi-allelic sites
  vcfTemp<-vcfTemp[vcfTemp$GT_TUMOR!="./." & vcfTemp$GT_TUMOR!="0/0",]
  vcfTemp <- vcfTemp[as.numeric(gsub("\\/.*", "", vcfTemp$GT_TUMOR)) == 0 |
                       as.numeric(gsub("\\/.*", "", vcfTemp$GT_TUMOR)) == as.numeric(gsub(".*\\/", "", vcfTemp$GT_TUMOR)) , ]
  
  # proper-view
  gtIndex <- as.numeric(gsub(".*\\/", "", vcfTemp$GT_TUMOR))
  
  toSplit <- c(
    "AF",
    "ALT",
    "AO",
    "FAO",
    "FDVR",
    "FR",
    "FSAF",
    "FSAR",
    "FWDB",
    "HRUN",
    "LEN",
    "MLLD",
    "OALT",
    "OID",
    "OMAPALT",
    "OPOS",
    "OREF",
    "PB",
    "PBP",
    "RBI",
    "REFB",
    "REVB",
    "SAF",
    "SAR",
    "SSEN",
    "SSEP",
    "SSSB",
    "STB",
    "STBP",
    "TYPE",
    "VARB",
    "FAO_NORMAL",
    "AF_NORMAL"
  )
  toSplit<-toSplit[toSplit%in%colnames(vcfTemp)]
  for (i in 1:length(toSplit)) {
    xTemp <- strsplit(vcfTemp[, colnames(vcfTemp) == toSplit[i]], ",")
    for (j in 1:length(xTemp)) {
      xTemp[[j]] <- xTemp[[j]][gtIndex[j]]
    }
    vcfTemp[, colnames(vcfTemp) == toSplit[i]] <- unlist(xTemp)
  }
  
  vcfTemp$MUTID<-paste(vcfTemp$CHROM,vcfTemp$POS,vcfTemp$REF,vcfTemp$ALT,sep=":")
  
  # attach TSV annotation -------------------------------------------------
  
  tsvTemp$alt<-gsub(".*\\/","",tsvTemp$genotype)
  tsvTemp$MUTID<-paste(tsvTemp$X..locus,tsvTemp$ref,tsvTemp$alt,sep=":")
  vcfTemp<-merge(vcfTemp,tsvTemp,by="MUTID",all.x=T)
  vcfTemp<-vcfTemp[,!colnames(vcfTemp)%in%c("MUTID","X..locus","type"
                                           ,"ref","length","genotype",
                                           "filter","normal_genotype",
                                           "pvalue","coverage",
                                           "allele_coverage","iscn",
                                           "confidence","precision",
                                           "Subtype","Call",
                                           "no_call_reason","normal_genotype.1",
                                           "normal_coverage","normal_allele_coverage",
                                           "allele_ratio","X._frequency","alt")]
    write.table(vcfTemp,file=paste(path.out,"/",sampleName,"_vcf_tsv_merged_properview.txt",sep=""),sep="\t",row.names=F,quote=F)
  } else {
    vcfFile<-grep(sampleName,list.files(path.in,pattern="Filtered",full.names=T,recursive=T),value=T)
    tsvFile<-grep(sampleName,list.files(path.in,pattern="\\.tsv",full.names=T,recursive=T),value=T)
    
    vcfTemp<-read.table(vcfFile,sep="\t",as.is=T,comment.char="#",header=F,stringsAsFactors=F)
    colnames(vcfTemp)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","TUMOR")
    vcfTemp<-vcfTemp[vcfTemp$ALT!="<CNV>",]
    
    tsvTemp<-read.table(tsvFile,sep="\t",header=T,comment.char="",skip=2,stringsAsFactors = F,quote="")
    tsvTemp<-tsvTemp[tsvTemp$type!="CNV",]
    
    if(nrow(vcfTemp)!=nrow(tsvTemp)){
      stop("Error: different number of variants in vcf and tsv files")
    }
    
    # splitting of INFO field -------------------------------------------------
    
    uniqueField<-unique(unlist(lapply(sapply(vcfTemp$INFO,strsplit,";"),function(x)gsub("=.*","",x))))
    infoTemp<-strsplit(vcfTemp$INFO,";")
    list2df<-function(x){
      y<-data.frame(field=gsub("=.*","",x),value=gsub(".*=","",x),stringsAsFactors = F)
    }
    infoTemp<-lapply(infoTemp,list2df)
    
    for(k in 1:length(infoTemp)){
      field=uniqueField[!uniqueField%in%infoTemp[[k]][,1]]
      if(length(field)>0){
        missing<-data.frame(field=field,value="")
        infoTemp[[k]]<-rbind(infoTemp[[k]],missing)
      }
      infoTemp[[k]]<-infoTemp[[k]][order(infoTemp[[k]][,1]),]
      infoTemp[[k]]<-infoTemp[[k]][infoTemp[[k]][,1]!="FUNC",]
    }
    
    if(length(table(sapply(infoTemp,nrow)))>1){
      stop("Error: the number of field in the INFO column is different across variants")
    } 
    columnNames<-as.vector(infoTemp[[1]]$field)
    infoTemp<-lapply(infoTemp,function(x)as.vector(x$value))
    dfInfo <- data.frame(matrix(unlist(infoTemp), nrow=length(infoTemp), byrow=T),stringsAsFactors = F)
    colnames(dfInfo)<-columnNames
    
    # attach to vcf data frame
    if(nrow(vcfTemp)!=nrow(dfInfo)){
      stop("Error: different number of rows between vcf and INFO data frame")
    }
    vcfTemp<-cbind(vcfTemp,dfInfo)
    
    # splitting of TUMOR field -------------------------------------------------
    
    tumorTemp<-lapply(strsplit(vcfTemp$TUMOR,":"),"[",1:2)
    dfTumor<-data.frame(GT_TUMOR=sapply(tumorTemp,"[",1),GQ_TUMOR=sapply(tumorTemp,"[",2),stringsAsFactors = F)
    
    # splitting of NORMAL field -------------------------------------------------
    
    # attach to vcf data frame
    if(nrow(vcfTemp)!=nrow(dfTumor)){
      stop("Error: different number of rows between vcf and TUMOR data frame")
    }
    vcfTemp<-cbind(vcfTemp,dfTumor)
    
    vcfTemp<-vcfTemp[,!colnames(vcfTemp)%in%c("INFO","FORMAT","TUMOR")]
	
	# remove multi-allelic sites
  vcfTemp<-vcfTemp[vcfTemp$GT_TUMOR!="./." & vcfTemp$GT_TUMOR!="0/0",]
  vcfTemp <- vcfTemp[as.numeric(gsub("\\/.*", "", vcfTemp$GT_TUMOR)) == 0 |
                       as.numeric(gsub("\\/.*", "", vcfTemp$GT_TUMOR)) == as.numeric(gsub(".*\\/", "", vcfTemp$GT_TUMOR)) , ]
  
  # proper-view
  gtIndex <- as.numeric(gsub(".*\\/", "", vcfTemp$GT_TUMOR))
  
  toSplit <- c(
    "AF",
    "ALT",
    "AO",
    "FAO",
    "FDVR",
    "FR",
    "FSAF",
    "FSAR",
    "FWDB",
    "HRUN",
    "LEN",
    "MLLD",
    "OALT",
    "OID",
    "OMAPALT",
    "OPOS",
    "OREF",
    "PB",
    "PBP",
    "RBI",
    "REFB",
    "REVB",
    "SAF",
    "SAR",
    "SSEN",
    "SSEP",
    "SSSB",
    "STB",
    "STBP",
    "TYPE",
    "VARB",
    "FAO_NORMAL",
    "AF_NORMAL"
  )
  toSplit<-toSplit[toSplit%in%colnames(vcfTemp)]
  for (i in 1:length(toSplit)) {
    xTemp <- strsplit(vcfTemp[, colnames(vcfTemp) == toSplit[i]], ",")
    for (j in 1:length(xTemp)) {
      xTemp[[j]] <- xTemp[[j]][gtIndex[j]]
    }
    vcfTemp[, colnames(vcfTemp) == toSplit[i]] <- unlist(xTemp)
  }
	
    vcfTemp$MUTID<-paste(vcfTemp$CHROM,vcfTemp$POS,vcfTemp$REF,vcfTemp$ALT,sep=":")
    
    # attach TSV annotation -------------------------------------------------
    
    tsvTemp$alt<-gsub(".*\\/","",tsvTemp$genotype)
    tsvTemp$MUTID<-paste(tsvTemp$X..locus,tsvTemp$ref,tsvTemp$alt,sep=":")
    vcfTemp<-merge(vcfTemp,tsvTemp,by="MUTID",all.x=T)
    vcfTemp<-vcfTemp[,!colnames(vcfTemp)%in%c("MUTID","X..locus","type"
                                              ,"ref","length","genotype",
                                              "filter",
                                              "pvalue","coverage",
                                              "allele_coverage","alt")]
    write.table(vcfTemp,file=paste(path.out,"/",sampleName,"_vcf_tsv_merged_properview.txt",sep=""),sep="\t",row.names=F,quote=F)
  }
    }
  




