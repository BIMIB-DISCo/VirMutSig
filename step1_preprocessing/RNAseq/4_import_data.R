# load reference data
load("reference.RData")

# list files
vcf_files = list.files("vcfDir",recursive=TRUE)
coverage_files = list.files("coverage",recursive=TRUE)

# read variants
variants = NULL
for(f in paste0("vcfDir/",sort(vcf_files))) {
    if(file.info(f)$size>0) {
        curr_data = read.table(file=f,header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
        curr_data = curr_data[,c("Position","Ref","VarAllele","VarAllele","Reads2","Reads1","VarFreq","Pvalue")]
        curr_data$Reads1 = curr_data$Reads1+curr_data$Reads2
        curr_patient = rep(strsplit(f,split="/")[[1]][2],nrow(curr_data))
        curr_variants = paste0(curr_data[,1],"_",curr_data[,2],"_",curr_data[,3])
        curr_data = cbind(curr_patient,curr_variants,curr_data)
        colnames(curr_data) = c("SampleId","VariantId","Position","Reference","Variant","ReferenceTrinucleotide","VariantCount","TotalCount","VariantFrequency","Pvalue")
        curr_data$VariantFrequency = (as.numeric(gsub("%","",curr_data$VariantFrequency))/100)
        curr_data$ReferenceTrinucleotide = paste0(reference[(curr_data$Position-1)],curr_data$ReferenceTrinucleotide,reference[(curr_data$Position+1)])
        variants = rbind(variants,curr_data)
    }
}

# read coverage files
metadata = array(NA,c(length(coverage_files),4))
colnames(metadata) = c("MeanCoverage","MedianCoverage","PercentageHighCoverage","TotalCoverage")
cont = 0
patients = NULL
for(f in paste0("coverage/",sort(coverage_files))) {
    if(file.info(f)$size>0) {
        curr_data = read.table(file=f,header=FALSE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
        patients = c(patients,gsub(".txt","",strsplit(f,split="/")[[1]][[2]]))
        cont = cont + 1
        metadata[cont,"MeanCoverage"] = mean(curr_data[,3])
        metadata[cont,"MedianCoverage"] = median(curr_data[,3])
        metadata[cont,"PercentageHighCoverage"] = length(which(curr_data[,3]>=50))/nrow(curr_data)
        metadata[cont,"TotalCoverage"] = sum(as.numeric(curr_data[,3]))
    }
}
metadata = metadata[1:length(patients),]
rownames(metadata) = patients

# save results
save(variants,file="results/variants.RData")
save(metadata,file="results/metadata.RData")
