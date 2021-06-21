library("optparse")

option_list = list(
    make_option(c("-m", "--SNV_list"), type="character", default="SNV_list.txt", 
                help="Input of SNV list obtained after preprocessing [default=%default]", metavar="character"),
    make_option(c("-c", "--config_file"), type="character", default="denovo_config.yaml", 
                help="specify the YAML config file [default=%default]", metavar="character"),
    make_option(c("-o", "--output_dir"), type="character", default="denovo_results", 
                help="output directory name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# install/load required libraries and scripts

# Package names
packages <- c("parallel","ggplot2", "data.table", "gridExtra", "NMF", "nnls", "yaml")

# rimuoverlo e metterlo nel readme su github
# Install packages not yet installed
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#     install.packages(packages[!installed_packages])
# }

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

parameters = yaml::read_yaml(opt$config_file)

source("denovo_utils.R")

variants = read.table(opt$SNV_list,header=TRUE,sep=";",check.names=FALSE,stringsAsFactors=FALSE)

# preprocessing
invalid_samples = which(table(variants$SampleId)>as.numeric(parameters$MAX_SNV_SAMPLE))
if(length(invalid_samples)>0) {
    invalid_samples = names(invalid_samples)
    variants = variants[which(!variants$SampleId%in%invalid_samples),]
}
if(parameters$MINOR_SNV_SEL=="always") { # consider always minor variants
    always_minor = sort(unique(variants$VariantId))
    always_minor_freq = rep(NA,length(always_minor))
    cont = 0
    for(i in always_minor) {
        cont = cont + 1
        always_minor_freq[cont] = max(variants$VariantFrequency[which(variants$VariantId==i)])
    }
    valid_variants = always_minor[which(always_minor_freq<parameters$CLONAL_SNV_THR)]
    variants = variants[which(variants$VariantId%in%valid_variants),]
}
if(parameters$MINOR_SNV_SEL=="all") { # consider all minor variants
    variants = variants[which(variants$VariantFrequency<parameters$CLONAL_SNV_THR),]
}
invalid_samples = which(table(variants$SampleId)<as.numeric(parameters$MIN_SNV_SAMPLE))
if(length(invalid_samples)>0) {
    invalid_samples = names(invalid_samples)
    variants = variants[which(!variants$SampleId%in%invalid_samples),]
}

# build variants matrix
if(as.numeric(parameters$NUM_CONTEXTS)==6) { # 6 contexts analysis
    variants_matrix = array(0,c(length(unique(variants$SampleId)),6))
    rownames(variants_matrix) = sort(unique(variants$SampleId))
    colnames(variants_matrix) = c("C>A:G>T","C>G:G>C","C>T:G>A","T>A:A>T","T>C:A>G","T>G:A>C")
    for(i in colnames(variants_matrix)) {
        curr_variants = strsplit(strsplit(i,":")[[1]],">")
        curr_variants1 = curr_variants[[1]]
        curr_variants2 = curr_variants[[2]]
        curr_variants = sort(c(which(variants$Reference==curr_variants1[1]&variants$Alternative==curr_variants1[2]),
                               which(variants$Reference==curr_variants2[1]&variants$Alternative==curr_variants2[2])))
        curr_variants = table(variants$SampleId[curr_variants])
        variants_matrix[names(curr_variants),i] = as.numeric(curr_variants)
    }
}

if(as.numeric(parameters$NUM_CONTEXTS)==96) { # 96 contexts analysis
    variants_matrix = array(0,c(length(unique(variants$SampleId)),96))
    rownames(variants_matrix) = sort(unique(variants$SampleId))
    colnames(variants_matrix) = trinucleotides
    variants_trinucleotides = variants$ReferenceTrinucleotide
    for(i in 1:length(variants_trinucleotides)) {
        variants_trinucleotides[i] = paste0(substr(variants_trinucleotides[i],1,1),
                                            "[",paste0(variants[i,"Reference"],
                                                       ">",
                                                       variants[i,"Variant"]),
                                            "]",
                                            substr(variants_trinucleotides[i],3,3))
    }
    variants_trinucleotides = gsub("G>T","C>A",variants_trinucleotides)
    variants_trinucleotides = gsub("G>C","C>G",variants_trinucleotides)
    variants_trinucleotides = gsub("G>A","C>T",variants_trinucleotides)
    variants_trinucleotides = gsub("A>T","T>A",variants_trinucleotides)
    variants_trinucleotides = gsub("A>G","T>C",variants_trinucleotides)
    variants_trinucleotides = gsub("A>C","T>G",variants_trinucleotides)
    for(i in colnames(variants_matrix)) {
        curr_variants = which(variants_trinucleotides==i)
        if(length(curr_variants)>0) {
            curr_variants = table(variants$SampleId[curr_variants])
            variants_matrix[names(curr_variants),i] = as.numeric(curr_variants)
        }
    }
}

# perform inference
num_processes = as.numeric(parameters$N_CORE)
if(num_processes==0) {
    num_processes = Inf
}
seed = as.numeric(parameters$SEED)
if(seed==0) {
    seed = NULL
}
signatures_decomposition = signatures.decomposition(x=variants_matrix,K=as.numeric(parameters$MIN_NUM_SIG):as.numeric(parameters$MAX_NUM_SIG),
                                                    nmf_runs=as.numeric(parameters$NMF_ITER),
                                                    num_processes=num_processes,seed=seed)

# save results
dir.create(paste0(opt$output_dir,"/inputs"), recursive = T, showWarnings = F)

dir.create(paste0(opt$output_dir,"/rank/figures"), recursive = T, showWarnings = F)

dir.create(paste0(opt$output_dir,"/signatures/files"), recursive = T, showWarnings = F)
dir.create(paste0(opt$output_dir,"/signatures/figures"), recursive = T, showWarnings = F)

write.table(variants_matrix,file=paste0(opt$output_dir,"/inputs/sample_contest_count.txt"),
            quote=FALSE,
            sep=";",
            row.names=TRUE,
            col.names=TRUE)

# rank estimation
# write.table(as.numeric(signatures_decomposition$measures[,"Cophenetic_Coefficient"]),
#             file=paste0(opt$output_dir,"/rank/files/cophenetic_coefficient.txt"),
#             quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)
# 
# write.table(as.numeric(signatures_decomposition$measures[,"Dispersion_Coefficient"]),
#             file=paste0(opt$output_dir,"/rank/files/dispersion_coefficient.txt"),
#             quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)
# 
# write.table(as.numeric(signatures_decomposition$measures[,"Silhouette_Consensus"]),
#             file=paste0(opt$output_dir,"/rank/files/silhouette_consensus.txt"),
#             quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)
# 
# write.table(as.numeric(signatures_decomposition$measures[,"Explained_Variance"]),
#             file=paste0(opt$output_dir,"/rank/files/explained_variance.txt"),
#             quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)

write.table(x = signatures_decomposition$measures[, c("Rank","Cophenetic_Coefficient", "Dispersion_Coefficient", "Silhouette_Consensus", "Explained_Variance")], 
            file=paste0(opt$output_dir,"/rank/metric_coefficients.txt"), 
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE, 
            sep = ';')

RANGE = as.numeric(parameters$MIN_NUM_SIG):as.numeric(parameters$MAX_NUM_SIG)

VALUE = as.numeric(signatures_decomposition$measures[,"Silhouette_Consensus"])
p = ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + 
    geom_line(linetype="dashed") + 
    geom_point() + 
    scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + 
    ggtitle("Evaluation of Optimal Number of Signatures") + 
    xlab("Number of Signatures") + 
    ylab("Silhouette Consensus") + 
    theme_bw()
ggsave(paste0(opt$output_dir,"/rank/figures/silhouette_consensus.pdf"),plot=p,width=10,height=8)

VALUE = as.numeric(signatures_decomposition$measures[,"Cophenetic_Coefficient"])
p = ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + 
    geom_line(linetype="dashed") + 
    geom_point() + 
    scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + 
    ggtitle("Evaluation of Optimal Number of Signatures") + 
    xlab("Number of Signatures") +
    ylab("Cophenetic Coefficient") + 
    theme_bw()
ggsave(paste0(opt$output_dir,"/rank/figures/cophenetic_coefficient.pdf"),plot=p,width=10,height=8)

VALUE = as.numeric(signatures_decomposition$measures[,"Dispersion_Coefficient"])
p = ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + 
    geom_line(linetype="dashed") + 
    geom_point() + 
    scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + 
    ggtitle("Evaluation of Optimal Number of Signatures") + 
    xlab("Number of Signatures") + 
    ylab("Dispersion Coefficient") + 
    theme_bw()
ggsave(paste0(opt$output_dir,"/rank/figures/dispersion_coefficient.pdf"),plot=p,width=10,height=8)

VALUE = as.numeric(signatures_decomposition$measures[,"Explained_Variance"])
p = ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + 
    geom_line(linetype="dashed") + 
    geom_point() + 
    scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + 
    ggtitle("Evaluation of Optimal Number of Signatures") + 
    xlab("Number of Signatures") + 
    ylab("Explained Variance") + 
    theme_bw()
ggsave(paste0(opt$output_dir,"/rank/figures/explained_variance.pdf"),plot=p,width=10,height=8)

# de novo signatures
for(i in names(signatures_decomposition$alpha)) {
    curr_beta = signatures_decomposition$beta[[i]]
    
    print_beta <- as.data.frame(curr_beta)
    print_beta$signature <- rownames(print_beta)
    
    write.table(print_beta[,c(ncol(print_beta), 1:(ncol(print_beta)-1))],file=paste0(opt$output_dir,"/signatures/files/",i,".txt"),
                quote=FALSE,sep=";",row.names=FALSE,col.names=TRUE)
    
    if(as.numeric(parameters$NUM_CONTEXTS)==6) {
        p = plot.6.signatures(curr_beta)
        ggsave(paste0(opt$output_dir,"/signatures/figures/",i,".pdf"),
               plot=p,width=10,height=8)
        dev.off()
    }
    if(as.numeric(parameters$NUM_CONTEXTS)==96) {
        p = plot.96.signatures(curr_beta)
        ggsave(paste0(opt$output_dir,"/signatures/figures/",i,".pdf"),plot=p,width=15,height=12)
        dev.off()
        
    }
}
