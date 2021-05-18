# load required libraries and scripts
library("data.table")
library("ggplot2")
library("gridExtra")
library("NMF")
library("nnls")
source("R/signatures.analysis.R")

# function to plot 6-contexts signatures
"plot.6.signatures" = function(beta) {
    x <- as.data.table(melt(beta, varnames = c("signature", "cat")))
    glist <- list()
    for (i in 1:nrow(beta)) {
        plt <- ggplot(x[signature == rownames(beta)[i]]) + geom_bar(aes(x = cat, 
            y = value, fill = cat), stat = "identity", position = "identity") + 
            facet_wrap(~cat, nrow = 1, scales = "free_x") + theme(axis.text.x = element_text(angle = 90, 
            hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
            ggtitle(rownames(beta)[i]) + theme(text = element_text(size=12),legend.position = "none") + 
            ylab("Frequency of mutations") + xlab("") + ylim(0,1)
        plt <- plt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        glist[[i]] <- plt
    }
    p = grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))
    return(p)
}

# function to plot 96-contexts signatures
"plot.96.signatures" = function(beta) {
    x <- as.data.table(melt(beta, varnames = c("signature", "cat")))
    x[,Context:=paste0(substr(cat,1,1),".",substr(cat,7,7))]
    x[,alt:=paste0(substr(cat,3,3),">",substr(cat,5,5))]
    glist <- list()
    for (i in 1:nrow(beta)) {
        plt <- ggplot(x[signature==rownames(beta)[i]]) + 
                geom_bar(aes(x=Context,y=value,fill=alt),stat="identity",position="identity") + 
                facet_wrap(~alt,nrow=1,scales="free_x") + 
                theme(axis.text.x=element_text(angle=90,hjust=1),panel.background=element_blank(),axis.line=element_line(colour="black")) + 
                ggtitle(rownames(beta)[i]) + theme(legend.position="none") + ylab("Frequency of mutations")
        plt <- plt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        glist[[i]] <- plt
    }
    p = grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))
    return(p)
}

# read data
load(file="trinucleotides.RData")
variants = read.table("SNV_list.txt",header=TRUE,sep=";",check.names=FALSE,stringsAsFactors=FALSE)
denovo_config = read.table("denovo_config.txt",header=FALSE,sep=";",check.names=FALSE,stringsAsFactors=FALSE,row.names=1)

# preprocessing
invalid_samples = which(table(variants$SampleId)>as.numeric(denovo_config["MAX_SNV_SAMPLE",]))
if(length(invalid_samples)>0) {
    invalid_samples = names(invalid_samples)
    variants = variants[which(!variants$SampleId%in%invalid_samples),]
}
if(denovo_config["MINOR_SNV_SEL",]=="always") { # consider always minor variants
    always_minor = sort(unique(variants$VariantId))
    always_minor_freq = rep(NA,length(always_minor))
    cont = 0
    for(i in always_minor) {
        cont = cont + 1
        always_minor_freq[cont] = max(variants$VariantFrequency[which(variants$VariantId==i)])
    }
    valid_variants = always_minor[which(always_minor_freq<denovo_config["CLONAL_SNV_THR",])]
    variants = variants[which(variants$VariantId%in%valid_variants),]
}
if(denovo_config["MINOR_SNV_SEL",]=="all") { # consider all minor variants
    variants = variants[which(variants$VariantFrequency<denovo_config["CLONAL_SNV_THR",]),]
}
invalid_samples = which(table(variants$SampleId)<as.numeric(denovo_config["MIN_SNV_SAMPLE",]))
if(length(invalid_samples)>0) {
    invalid_samples = names(invalid_samples)
    variants = variants[which(!variants$SampleId%in%invalid_samples),]
}

# build variants matrix
if(as.numeric(denovo_config["NUM_CONTEXTS",])==6) { # 6 contexts analysis
    variants_matrix = array(0,c(length(unique(variants$SampleId)),6))
    rownames(variants_matrix) = sort(unique(variants$SampleId))
    colnames(variants_matrix) = c("C>A:G>T","C>G:G>C","C>T:G>A","T>A:A>T","T>C:A>G","T>G:A>C")
    for(i in colnames(variants_matrix)) {
        curr_variants = strsplit(strsplit(i,":")[[1]],">")
        curr_variants1 = curr_variants[[1]]
        curr_variants2 = curr_variants[[2]]
        curr_variants = sort(c(which(variants$Reference==curr_variants1[1]&variants$Variant==curr_variants1[2]),which(variants$Reference==curr_variants2[1]&variants$Variant==curr_variants2[2])))
        curr_variants = table(variants$SampleId[curr_variants])
        variants_matrix[names(curr_variants),i] = as.numeric(curr_variants)
    }
}
if(as.numeric(denovo_config["NUM_CONTEXTS",])==96) { # 96 contexts analysis
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
num_processes = as.numeric(denovo_config["N_CORE",])
if(num_processes==0) {
    num_processes = Inf
}
seed = as.numeric(denovo_config["SEED",])
if(seed==0) {
    seed = NULL
}
signatures_decomposition = signatures.decomposition(x=variants_matrix,K=as.numeric(denovo_config["MIN_NUM_SIG",]):as.numeric(denovo_config["MAX_NUM_SIG",]),nmf_runs=as.numeric(denovo_config["NMF_ITER",]),num_processes=num_processes,seed=seed)

# save results
dir.create("results")
dir.create("results/inputs")
dir.create("results/rank")
dir.create("results/rank/files")
dir.create("results/rank/figures")
dir.create("results/signatures")
dir.create("results/signatures/files")
dir.create("results/signatures/figures")
write.table(variants_matrix,file="results/inputs/variants_matrix.txt",quote=FALSE,sep=";",row.names=TRUE,col.names=TRUE)

# rank estimation
write.table(as.numeric(signatures_decomposition$measures[,"Cophenetic_Coefficient"]),file="results/rank/files/cophenetic_coefficient.txt",quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)
write.table(as.numeric(signatures_decomposition$measures[,"Dispersion_Coefficient"]),file="results/rank/files/dispersion_coefficient.txt",quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)
write.table(as.numeric(signatures_decomposition$measures[,"Silhouette_Consensus"]),file="results/rank/files/silhouette_consensus.txt",quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)
write.table(as.numeric(signatures_decomposition$measures[,"Explained_Variance"]),file="results/rank/files/explained_variance.txt",quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)
RANGE = as.numeric(denovo_config["MIN_NUM_SIG",]):as.numeric(denovo_config["MAX_NUM_SIG",])
VALUE = as.numeric(signatures_decomposition$measures[,"Silhouette_Consensus"])
p = ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Silhouette Consensus")
ggsave("results/rank/figures/silhouette_consensus.pdf",plot=p,width=20,height=15)
VALUE = as.numeric(signatures_decomposition$measures[,"Cophenetic_Coefficient"])
p = ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Cophenetic Coefficient")
ggsave("results/rank/figures/cophenetic_coefficient.pdf",plot=p,width=20,height=15)
VALUE = as.numeric(signatures_decomposition$measures[,"Dispersion_Coefficient"])
p = ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Dispersion Coefficient")
ggsave("results/rank/figures/dispersion_coefficient.pdf",plot=p,width=20,height=15)
VALUE = as.numeric(signatures_decomposition$measures[,"Explained_Variance"])
p = ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Explained Variance")
ggsave("results/rank/figures/explained_variance.pdf",plot=p,width=20,height=15)

# de novo signatures
for(i in names(signatures_decomposition$alpha)) {
    curr_beta = signatures_decomposition$beta[[i]]
    write.table(curr_beta,file=paste0("results/signatures/files/",i,".txt"),quote=FALSE,sep=";",row.names=TRUE,col.names=FALSE)
    if(as.numeric(denovo_config["NUM_CONTEXTS",])==6) {
        p = plot.6.signatures(curr_beta)
        ggsave(paste0("results/signatures/figures/",i,".pdf"),plot=p,width=20,height=15)
        dev.off()
    }
    if(as.numeric(denovo_config["NUM_CONTEXTS",])==96) {
        p = plot.96.signatures(curr_beta)
        ggsave(paste0("results/signatures/figures/",i,".pdf"),plot=p,width=20,height=15)
        dev.off()

    }
}
