suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(fitdistrplus))

option_list = list(
  make_option(
    "--seuratObj",
    dest="seurat_file_path",
    default=NA,
    type='character',
    help="File path to the rds file containing the seurat object."),
  make_option("--output",
              default=NA,
              type='character',
              help="Output name"),
  make_option(
	"--scriptFunctions",
	default=NA,
	type='character',
	help="path of the script_funtions_covid.R file"
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

source(opt$scriptFunctions)

seuratObj <- readRDS(opt$seurat_file_path)
#############################
########## METADATA
#############################
set.seed(20200221)

print("########################## HTO metaData ##########################")
# creating of HTO.data, the metadata
HTO.data <- as.data.frame(Matrix::t(seuratObj@assays$HTO@counts)) #hashcounts
HTO.data.raw <- HTO.data
HTO.data$BARCODE <- rownames(HTO.data)
features <- colnames(HTO.data.raw) #hashtags
HTO.data$maxHTO <- apply(HTO.data.raw,1,which.max) #max hash
HTO.data <- mutate(HTO.data,maxHTO=features[maxHTO])
HTO.data$maxval <- apply(HTO.data.raw,1,max)
sec.fun <- function(x){features[order(x)[length(x)-1]]}
secval.fun <- function(x){x[order(x)[length(x)-1]]}
HTO.data$secHTO <- apply(HTO.data.raw,1,sec.fun)
HTO.data$secval <- apply(HTO.data.raw,1,secval.fun)
HTO.data$nUMI <- seuratObj$nCount_RNA
HTO.data$nGene <- seuratObj$nFeature_RNA
HTO.data$cell <- rownames(HTO.data.raw)
# add new info to the metadata: w.crit / dif.crit / mix.crit
HTO.data <- dplyr::group_by(HTO.data,maxHTO) %>%
  dplyr::mutate(w.crit=within.crit(maxval,q=0.999)) %>%
  dplyr::mutate(dif.crit=dif.crit(maxval,secval,p=0.8)) %>%
  dplyr::mutate(mix.crit=mix.crit(secval,q=0.99)) %>%
  dplyr::ungroup()

# creating the corresponding labels (for singlets what is the label)
HTO.data$dif_label <- ifelse(HTO.data$dif.crit=="ok",HTO.data$maxHTO,HTO.data$dif.crit)
HTO.data$mix_label <- ifelse(HTO.data$mix.crit=="Singlet",HTO.data$maxHTO,HTO.data$mix.crit)
HTO.data$comp_label <-ifelse(HTO.data$w.crit=="Negative","Negative",
                             ifelse(HTO.data$dif.crit=="too.close","Doublet",HTO.data$mix_label))

seuratObj$HTOmaxval <- HTO.data$maxval
seuratObj@misc$HTO.data <- HTO.data

if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,"_SEURAT__HTO_METADATA.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output))
}
