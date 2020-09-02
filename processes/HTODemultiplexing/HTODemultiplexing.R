suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(openxlsx))

option_list = list(
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="File path to the rds file containing the seurat object."
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="Output name"
  ),
  make_option(
    "--freemuxFile",
    default=NA,
    type='character',
    help="File path to the freemuxlet file."
  ),
  make_option(
    "--qcFile",
    default=NA,
    type='character',
    help="File path to the QC log file."
	),
	make_option(
      "--scriptFunctions",
      default=NA,
      type='character',
      help="path of the script_funtions_covid.R file"
  	)

)
opt <- parse_args(OptionParser(option_list=option_list))

source(opt$scriptFunctions)

dir.create("Robjects")
seuratObj <- readRDS(opt$seuratObj)

#FREEMUXLETDATA
freemuxlet <- fread(paste0("zcat < ", opt$freemuxFile),
                    data.table=F)
rownames(freemuxlet) <- freemuxlet$BARCODE

##### HTO DEMUX
########################################
# The demultiplexing function HTODemux() implements the following procedure:
# -We perform a k-medoid clustering on the normalized HTO values, which initially
#  separates cells into K(# of samples)+1 clusters.)
# -We calculate a ‘negative’ distribution for HTO. For each HTO, w
#  we use the cluster with the lowest average value as the negative group.
# -For each HTO, we fit a negative binomial distribution to the negative cluster.
#  We use the 0.99 quantile of this distribution as a threshold.
#  Based on these thresholds, each cell is classified as positive or negative for each HTO.
#  Cells that are positive for more than one HTOs are annotated as doublets.
#  Most important parameter for tuning is positive.quantile (0.99),
#  to be positive it should be higher than that
set.seed(20200221)

seuratObj <- subset(seuratObj, subset = HTOmaxval > 0)
HTO.data <- seuratObj@misc$HTO.data[seuratObj@misc$HTO.data$maxval >0,]

seuratObj <- HTODemux(seuratObj,
                      assay = "HTO",
                      positive.quantile = 0.99,
                      nsamples = length(grep("Hashtag",colnames(HTO.data))),
                      init = length(grep("Hashtag",colnames(HTO.data)))+1)

# add  it on the meta-data
HTO.data$Demux_maxID <- seuratObj$HTO_maxID
HTO.data$Demux_secondID <- seuratObj$HTO_secondID
HTO.data$Demux_margin <- seuratObj$HTO_margin
HTO.data$Demux_glob <- seuratObj$HTO_classification
HTO.data$Demux_class <- seuratObj$HTO_classification.global
HTO.data$Demux_class <- factor(HTO.data$Demux_class,
                               levels = c("Singlet","Doublet","Negative"),
                               labels = c("Singlet","Doublet","Negative"))
HTO.data$Demux_label <- seuratObj$hash.ID


##### MULTISEQ DEMUX
########################################
# They use the trimmed PDF (0.1-0.9%) and take mode for negative and highest local max as positive,
# next derive iterative thresholds that maximize number of singlets between these 2 peaks
# tunable parameters? You can do autoThresh instead of fixed default quantile of 0.7
# calculation
set.seed(20200221)
seuratObj <- tryCatch(MULTIseqDemux(seuratObj,
                                    assay = "HTO",
                                    autoThresh = T,
                                    maxiter=5,
                                    qrange = seq(from = 0.1, to =  0.9, by =0.05)),
                      error = function(e) {print(paste("MY_ERROR:  ",e))
                        seuratObj$MULTI_ID <- rep("Failed",nrow(seuratObj@meta.data))
                        seuratObj$MULTI_classification <- rep("Failed",nrow(seuratObj@meta.data))
                        seuratObj$Multi_label <- rep("Failed",nrow(seuratObj@meta.data))
                        return(seuratObj)
                      })

# adding it on the meta-data
HTO.data$Multi_label <- seuratObj$MULTI_ID
HTO.data$Multi_glob <- seuratObj$MULTI_classification
HTO.data$Multi_class <- HTO.data$Multi_label
levels(HTO.data$Multi_class) <- ifelse(levels(HTO.data$Multi_class) %in% c("Doublet","Negative","Failed"),
                                       levels(seuratObj$MULTI_ID),
                                       "Singlet")
if(length(levels(HTO.data$Multi_class))==3){
  HTO.data$Multi_class <- factor(HTO.data$Multi_class,levels(HTO.data$Multi_class)[c(2,1,3)])
}

HTO.data$w.crit<-factor(HTO.data$w.crit, levels = c("Singlet", "Doublet", "Negative"))
HTO.data$mix.crit<-factor(HTO.data$mix.crit, levels = c("Singlet", "Doublet"))

print("FREEMUXLET")
########################################
seuratObj@meta.data$freemuxlet_type <- freemuxlet[colnames(seuratObj), 'DROPLET.TYPE']
seuratObj@meta.data$freemuxlet_guess <- freemuxlet[colnames(seuratObj), 'BEST.GUESS']

freemuxlet <- freemuxlet[colnames(seuratObj),]

Freemuxlet_type <- mapvalues(as.factor(freemuxlet[colnames(seuratObj),'DROPLET.TYPE']),
                             from = c("AMB","DBL","SNG"),
                             to = c("AMB","Doublet","Singlet"))

Seurat_global_class <- seuratObj@meta.data[rownames(freemuxlet),which(colnames(seuratObj@meta.data)=="HTO_classification.global")]

#levels(Freemuxlet_type) <- c("AMB","Doublet","Singlet","Negative")
#levels(Seurat_global_class) <- c("AMB","Doublet","Singlet","Negative")

Freemuxlet_type <- base::factor(Freemuxlet_type,levels=c("AMB","Doublet","Singlet","Negative"))
Seurat_global_class <- base::factor(Seurat_global_class,levels=c("AMB","Doublet","Singlet","Negative"))

confusionMatrix_output <- confusionMatrix(Freemuxlet_type,Seurat_global_class)
TheMatrix <- data.frame(confusionMatrix_output$table)

coresTable <- table(seuratObj$HTO_maxID, seuratObj$freemuxlet_guess)
coresTable <- coresTable[,unlist(lapply(colnames(coresTable),function(x) strsplit(x,",")[[1]][1] == strsplit(x,",")[[1]][2]))]

rownames(coresTable) <- gsub("Hashtag","#",rownames(coresTable))
secLabel <- as.data.frame(apply(coresTable,1,function(x) colnames(coresTable)[which.max(x)]))

freemuxlet_label <- apply(seuratObj@meta.data,
                          1,
                          function(x) ifelse(x[colnames(seuratObj@meta.data) == "freemuxlet_type"] =="SNG",
                                                                   rownames(secLabel)[which(secLabel[,1] == x[colnames(seuratObj@meta.data) == "freemuxlet_guess"])],
                                                                   x[colnames(seuratObj@meta.data) == "freemuxlet_type"]))
# add it to the seuratobj
seuratObj@meta.data$freemuxlet_label <- freemuxlet_label
# add it to the HTO data
HTO.data$freemuxlet_label <- freemuxlet_label
seuratObj@misc$HTO.data <- HTO.data

seuratObj$subsample <-
  paste0(seuratObj@project.name,".",
         ifelse(HTO.data$Demux_label=="Negative" & HTO.data$Multi_label=="Negative",
                "Negative",
                ifelse(HTO.data$Demux_label=="Doublet" & HTO.data$Multi_label=="Doublet",
                       "DBL",
                       freemuxlet_label)))

seuratObj@misc$TheMatrix <- TheMatrix

# OUTPUTS
write.xlsx(confusionMatrix_output$table, file = paste0(seuratObj@project.name,"_confusionmatrix.xlsx"))

qcfilename <- ifelse(is.na(opt$qcFile), paste0(seuratObj@project.name,"_logQC.txt"), opt$qcFile)
write("#QC: ~~~~DEMULTIPLEXING~~~~ < Intermediate QC",file = qcfilename,append=T)
write.table(head(seuratObj@meta.data,20),file = qcfilename, append = T,row.names = F)
write("#QC: ~~~~DEMULTIPLEXING~~~~ > Intermediate QC",file = qcfilename,append=T)

if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,"_SEURAT__HTO_DEMULTIPLEXING.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output))
}

