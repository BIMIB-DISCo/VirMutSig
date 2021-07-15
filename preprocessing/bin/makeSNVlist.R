#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

SRRfiles <- strsplit(x=args[1],split = " ",fixed=TRUE)[[1]]
DEPTHfiles <- strsplit(x=args[2],split = " ",fixed=TRUE)[[1]]

refFASTAf <- file(args[3], "r")

pv_thr <- as.numeric(strsplit(x=args[4], split=':', fixed=TRUE)[[1]][2])
var_freq_thr <- as.numeric(strsplit(x=args[5], split=':', fixed=TRUE)[[1]][2])
min_cov <- as.numeric(strsplit(x=args[6], split=':', fixed=TRUE)[[1]][2])
alt_read_thr <- as.numeric(strsplit(x=args[7], split=':', fixed=TRUE)[[1]][2])

refFASTA <- c()

for (l in readLines(con=refFASTAf,warn=FALSE)) {
  
  if(!startsWith(l,">")) {
    l <- gsub("[[:space:]]", "", l)
    nucleotides <- strsplit(x=l,split="",fixed=TRUE)[[1]]
    refFASTA <- c(refFASTA,nucleotides) 
  }  
}

# read variants
variants = NULL
for(SRRf in SRRfiles) {
    if(file.info(SRRf)$size>0) {
        curr_data = read.table(file=SRRf,header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
        curr_data = curr_data[,c("Position","Ref","VarAllele","Ref","Reads2","Reads1","VarFreq","Pvalue")]
        curr_data$Reads1 = curr_data$Reads1+curr_data$Reads2
        curr_patient = rep(strsplit(SRRf,split="\\.")[[1]][1],nrow(curr_data))
        curr_variants = paste0(curr_data[,1],"_",curr_data[,2],"_",curr_data[,3])
        curr_data = cbind(curr_patient,curr_variants,curr_data)
        colnames(curr_data) = c("SampleId","VariantId","Position","Reference","Alternative","ReferenceTrinucleotide","VariantCount","TotalCount","VariantFrequency","Pvalue")
        
        curr_data$VariantFrequency = (as.numeric(gsub("%","",curr_data$VariantFrequency))/100)
        curr_data$ReferenceTrinucleotide = paste0(refFASTA[(curr_data$Position-1)],curr_data$ReferenceTrinucleotide,refFASTA[(curr_data$Position+1)])
        variants = rbind(variants,curr_data)
    }
}

# read coverage files
metadata = array(NA,c(length(DEPTHfiles),4))
colnames(metadata) = c("MeanCoverage","MedianCoverage","PercentageHighCoverage","TotalCoverage")

cont = 0
patients = NULL
for(DEPTHf in DEPTHfiles) {
    if(file.info(DEPTHf)$size>0) {
        curr_data = read.table(file=DEPTHf,header=FALSE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
        patients = c(patients,strsplit(x=DEPTHf,split="\\.")[[1]][1])
        cont = cont + 1
        metadata[cont,"MeanCoverage"] = mean(curr_data[,3])
        metadata[cont,"MedianCoverage"] = median(curr_data[,3])
        metadata[cont,"PercentageHighCoverage"] = length(which(curr_data[,3]>=50))/nrow(curr_data)
        metadata[cont,"TotalCoverage"] = sum(as.numeric(curr_data[,3]))
    }
}


#saveRDS(variants, "variants.rds")

#metadata = metadata[1:length(patients),]
#rownames(metadata) = patients
#saveRDS(metadata, "metadata.rds")

# select only point mutations
processed_variants = variants[which(variants$Reference%in%c("A","C","G","T")&variants$Alternative%in%c("A","C","G","T")),]

# filering
idx_pv_thr = processed_variants$Pvalue<pv_thr
idx_var_freq_thr = processed_variants$VariantFrequency>=var_freq_thr
idx_min_cov = processed_variants$TotalCount>=min_cov
idx_alt_read_thr = processed_variants$VariantCount>=alt_read_thr
selected_processed_variants = processed_variants[which(idx_pv_thr & idx_var_freq_thr & idx_min_cov & idx_alt_read_thr),]

# write results to file
write.table(selected_processed_variants,file="SNV_list.txt",quote=FALSE,sep=";",row.names=FALSE,col.names=TRUE)
