# load data
load(file="results/variants.RData")
filtering_params = read.table("filtering_params.txt",sep=";",header=FALSE,row.names=1)

# select only point mutations
processed_variants = variants[which(variants$Reference%in%c("A","C","G","T")&variants$Variant%in%c("A","C","G","T")),]

# filering
pv_thr = variants$Pvalue<filtering_params["pv_thr",]
var_freq_thr = variants$VariantFrequency>=filtering_params["var_freq_thr",]
min_cov = variants$TotalCount>=filtering_params["min_cov",]
alt_read_thr = variants$VariantCount>=filtering_params["alt_read_thr",]
processed_variants = processed_variants[which(pv_thr&var_freq_thr&min_cov&alt_read_thr),]

# write results to file
write.table(processed_variants,file="SNV_list.txt",quote=FALSE,sep=";",row.names=FALSE,col.names=TRUE)
