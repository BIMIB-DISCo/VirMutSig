library("optparse")

option_list = list(
    make_option(c("-m", "--SNV_list"), type="character", default="SNV_list.txt", 
                help="Input of SNV list obtained after preprocessing [default=%default]", metavar="character"),
    make_option(c("-s", "--signatures"), type="character", default=NULL, 
                help="Input of chosen signature txt file (e.g., obtained after denovo analysis)", metavar="character"),
    make_option(c("-c", "--config_file"), type="character", default="assignment_config.yaml", 
                help="specify the YAML config file [default=%default]", metavar="character"),
    make_option(c("-o", "--output_file"), type="character", default="signature_assignment_result.txt", 
                help="output directory name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# load required libraries and scripts
packages <- c("parallel","NMF", "nnls", "lsa", "yaml")

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

source("assignment_utils.R")

variants = read.table(opt$SNV_list,header=TRUE,sep=";",check.names=FALSE,stringsAsFactors=FALSE)

signatures = read.table(opt$signatures,header=TRUE,sep=";",check.names=FALSE,stringsAsFactors=FALSE,row.names=1)

if(ncol(signatures)==6) {
    colnames(signatures) = c("C>A:G>T","C>G:G>C","C>T:G>A","T>A:A>T","T>C:A>G","T>G:A>C")
}
if(ncol(signatures)==96) {
    colnames(signatures) = trinucleotides
}

parameters = yaml::read_yaml(opt$config_file)

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
if(ncol(signatures)==6) { # 6 contexts analysis
    variants_matrix = array(0,c(length(unique(variants$SampleId)),6))
    rownames(variants_matrix) = sort(unique(variants$SampleId))
    colnames(variants_matrix) = c("C>A:G>T","C>G:G>C","C>T:G>A","T>A:A>T","T>C:A>G","T>G:A>C")
    for(i in colnames(variants_matrix)) {
        curr_variants = strsplit(strsplit(i,":")[[1]],">")
        curr_variants1 = curr_variants[[1]]
        curr_variants2 = curr_variants[[2]]
        curr_variants = sort(c(which(variants$Reference==curr_variants1[1]&variants$Alternative==curr_variants1[2]),which(variants$Reference==curr_variants2[1]&variants$Alternative==curr_variants2[2])))
        curr_variants = table(variants$SampleId[curr_variants])
        variants_matrix[names(curr_variants),i] = as.numeric(curr_variants)
    }
}
if(ncol(signatures)==96) { # 96 contexts analysis
    variants_matrix = array(0,c(length(unique(variants$SampleId)),96))
    rownames(variants_matrix) = sort(unique(variants$SampleId))
    colnames(variants_matrix) = trinucleotides
    variants_trinucleotides = variants$ReferenceTrinucleotide
    for(i in 1:length(variants_trinucleotides)) {
        variants_trinucleotides[i] = paste0(substr(variants_trinucleotides[i],1,1),"[",paste0(variants[i,"Reference"],">",variants[i,"Variant"]),"]",substr(variants_trinucleotides[i],3,3))
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
seed = as.numeric(parameters$SEED)
if(seed==0) {
    seed = NULL
}
set.seed(seed)
beta = signatures
alpha = array(NA,c(nrow(variants_matrix),nrow(signatures)))
rownames(alpha) = rownames(variants_matrix)
colnames(alpha) = rownames(signatures)
for(i in 1:nrow(alpha)) {
    alpha[i,] <- nnls(t(beta),as.vector(variants_matrix[i,]))$x
}
results = alpha
colnames(results) = paste0(colnames(results),"_exposure")
if(parameters$BOOTSTRAP=="yes") {
    N_CORE = as.numeric(parameters$N_CORE)
    if(N_CORE==0) {
        N_CORE = Inf
    }

    signatures_significance = signaturesSignificance(x=variants_matrix,
                                                     beta=beta,
                                                     cosine_thr=as.numeric(parameters$GOODNESS_FIT_THR),
                                                     min_contribution=as.numeric(parameters$MIN_SIG_FREQ),
                                                     pvalue_thr=as.numeric(parameters$PVALUE),
                                                     sparsify=FALSE,
                                                     nboot=as.numeric(parameters$NUM_ITER),
                                                     num_processes=N_CORE,
                                                     seed=seed,
                                                     verbose=TRUE)
    
    pvalues = signatures_significance$bootstrap_estimates$pvalues
    colnames(pvalues) = paste0(colnames(pvalues),"_pvalue")
    results = as.data.frame(cbind(results,pvalues))
}
results$sample <- rownames(results)
write.table(results[,c(ncol(results),1:(ncol(results)-1))],file=opt$output,quote=FALSE,sep=";",row.names=FALSE,col.names=TRUE)
