
parseIROutputOncomine<-function(file.path,
                                select.columns=c("sampleID","FUNC1.gene","INFO.A.TYPE","CHROM",
                                                 "INFO...OPOS","INFO...OREF","INFO...OALT","QUAL",
                                                 "INFO.A.MAF","INFO.1.DP","INFO.A.AO","INFO.1.MDP",
                                                 "INFO.A.MAO","FUNC1.function","FUNC1.protein",
                                                 "FUNC1.oncomineGeneClass","FUNC1.oncomineVariantClass",
                                                 "FUNC1.transcript","FUNC1.coding",
                                                 "FUNC1.location","FUNC1.polyphen","FUNC1.sift",
                                                 "dbsnp","maf","X5000Exomes","exac")){
  
  vcf.files<-list.files(file.path,recursive=T,"oncomine.tsv",full.names = T)
  annot.files<-list.files(file.path,recursive=T,"full.tsv",full.names = T)
  
  colonne<-vector("list",length=length(vcf.files))
  for(i in 1:length(vcf.files)){
    x<-try(read.table(vcf.files[i],header=T,sep="\t",as.is=T),silent=T)
    if(class(x)=="try-error"){
      print(paste("No mutations found for sample", gsub("_v1.*","",gsub(".*\\/","",vcf.files[i]))))
    } else {
      colonne[[i]]<-colnames(x)
    }
  }
  
  vcf.files<-vcf.files[sapply(colonne,length)>0]
  annot.files<-annot.files[sapply(colonne,length)>0]
  colonne<-colonne[sapply(colonne,length)>0]
  colonne<-Reduce(intersect,colonne)
  
  x<-read.table(vcf.files[1],header=T,sep="\t",as.is=T)
  x<-x[x$call=="POS",colnames(x)%in%colonne]
  x$sampleID<-gsub("_v1.*","",gsub(".*\\/","",vcf.files[1]))
  x$X..locus<-paste(x$CHROM,x$POS,sep=":")
  y<-read.table(annot.files[1],header=T,sep="\t",as.is=T,comment.char="",skip=2)
  if(identical(x$X..locus,y$X..locus)==FALSE){
    print("Mutations in vcf.file and annotation.file are differently ordered")
  }
  x<-cbind(x,y)
  for(i in 2:length(vcf.files)){
    temp<-read.table(vcf.files[i],header=T,sep="\t",as.is=T)
    temp<-temp[temp$call=="POS",colnames(temp)%in%colonne]
    temp$sampleID<-gsub("_v1.*","",gsub(".*\\/","",vcf.files[i]))
    temp$X..locus<-paste(temp$CHROM,temp$POS,sep=":")
    temp.annot<-read.table(annot.files[i],header=T,sep="\t",as.is=T,comment.char="",skip=2)
    if(identical(temp$X..locus,temp.annot$X..locus)==FALSE){
      print("Mutations in vcf.file and annotation.file are differently ordered")
    }
    temp<-cbind(temp,temp.annot)
    x<-rbind(x,temp)
  }
  x<-x[,select.columns]
  colnames(x)<-c("SampleID","Gene","Variant Type","Chr","Start","REF","ALT","QUAL",
                 "Molecular Allelic Frequency","Read Depth","Read Variant Count",
                 "Molecular Depth","Molecular Variant Count","Function","Protein Change",
                 "Oncomine Gene Class","Oncomine Variant Class","Transcript","Coding Change",
                 "Location","Plyphen","SIFT","dbSNP","MAF","5000 Exomes","ExAC")
  x$ExAC<-gsub("\\:.*","",gsub(".*AF_NFE=","",x$ExAC))
  x$ExAC[x$ExAC==""]<-NA
  return(x)
}
