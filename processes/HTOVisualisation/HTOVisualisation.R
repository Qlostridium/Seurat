suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggfan))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(patchwork))

option_list = list(
  make_option(
    "--seuratObj",
    dest="seurat_file_path",
    default=NA,
    type='character',
    help="File path to the rds file containing the seurat object."
  ),
  make_option(
    "--scriptFunctions",
    default=NA,
    type='character',
    help="path of the script_funtions_covid.R file"
  )
  
)
opt <- parse_args(OptionParser(option_list=option_list))
#opt <- list(seurat_file_path = "/home/quentinr/Documents/VIB_work/scRNA_seq/datasets/SAR13/work/72/5677d1deb7b2dd00ff761adcf6483c/SAR13.SEURAT__HTO_DEMULTIPLEXING.rds",
#            scriptFunctions = "/home/quentinr/Documents/VIB_work/scRNA_seq/Singularity_pipeline_submodules/Seurat/processes/utils/script_functions_COVID.R")
source(opt$scriptFunctions)

dir.create("Robjects")
dir.create("Plots/HTO",recursive = T)
seuratObj <- readRDS(opt$seurat_file_path)
#seuratObj <- readRDS("/home/quentinr/Documents/VIB_work/scRNA_seq/datasets/STA04/testRun/Robjects/seuratObj_demux_HTO.rds")
HTO.data <- seuratObj@misc$HTO.data
TheMatrix <- seuratObj@misc$TheMatrix
### METADATA
print("########################## visualisations ##########################")
# Cell frequencies per maximum Hashtag
max.barplot <-   ggplot(HTO.data,aes(x=maxHTO)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),fill = "#0073C2FF") +
  scale_y_continuous(labels=scales::percent) +
  labs(y = "Percent", fill="class") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=20)) +
  ggtitle("Cell frequencies per maximum Hashtag")
# Log2 counts per HTO
HTO_dens <- HTO.data[,grep("Hashtag",colnames(HTO.data))] %>%
  rownames_to_column(var="cell") %>%
  gather(key = HTO, value = counts,-cell) %>%
  ggplot(aes(x = log2(counts+1),color=HTO)) +
  geom_density() +
  theme_bw() +
  theme(text = element_text(size=15)) +
  facet_wrap(~HTO) +
  ggtitle("Log2 counts per HTO")
# Log2 counts per maximum and second largest value HTO
secmax_dens <- gather(HTO.data,c(maxval,secval), key = rank, value = counts) %>%
  ggplot(aes(x = log2(counts+1),color=rank)) +
  geom_density() +
  #xlim(0,count.lim) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  facet_wrap(~maxHTO) +
  ggtitle("Log2 counts per maximum and second largest value HTO")
secmax_mix <- ggplot(HTO.data,aes(x = log2(maxval+1),y = log2(secval+1),color=mix.crit)) +
  geom_point(size=0.75) +
  geom_abline(slope = 1 ,color = "blue" ) +
  geom_abline(slope = 0.80 ,color = "red" ) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  facet_wrap(~maxHTO) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(values=c("firebrick", "limegreen")) + #Doublet Singlet
  ggtitle("Log2 counts per maximum and second largest value HTO")
secmax_within <- ggplot(HTO.data,aes(x = log2(maxval+1),y = log2(secval+1),color=w.crit)) +
  geom_point(size=0.75) +
  geom_abline(slope = 1 ,color = "blue" ) +
  geom_abline(slope = 0.80 ,color = "red" ) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  facet_wrap(~maxHTO) +
  scale_color_manual(values=c("firebrick", "coral","limegreen")) + #Doublet Negative Singlet
  ggtitle("Log2 counts per maximum and second largest value HTO")
# MaxHTO
HTO_max <-dplyr::select(HTO.data,maxHTO,cell)
HTO_cnt <- HTO.data[,grep("Hashtag",colnames(HTO.data))] %>%
  rownames_to_column(var="cell") %>%
  gather(key = HTO, value = counts,-cell) %>%
  # mutate(HTO_num = as.numeric(as.factor(HTO))) %>%
  full_join(HTO_max,by="cell") %>%
  ggplot(aes(x=HTO,y=counts)) +
  geom_interval(aes(colour=..Interval..), intervals=c(0,0.5,0.75,0.9,0.99),size = 1.5) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  facet_wrap(~maxHTO,scales="free_y") +
  ggtitle("Distribution of counts split out by maxHTO") +
  xlab("Hashtag sample")
HTO_lines <- gather(HTO.data,grep("Hashtag",colnames(HTO.data),value = T), key = HTO, value = counts) %>%
  ggplot(aes(HTO, log2(counts+1), group=cell,color=maxHTO)) +
  geom_line(size=0.5,aes(color=maxHTO)) +
  geom_point() +
  theme(text = element_text(size=15)) +
  facet_wrap(~maxHTO,scales="free_y") +
  ggtitle("Line plot per maximum HTO")
#difference in raw counts
HTO_dif2 <- gather(HTO.data,c(maxval,secval), key = rank, value = counts) %>%
  ggplot(aes(rank, counts, group=cell,color=maxHTO)) +
  geom_line(size=0.5,aes(color=maxHTO)) +
  geom_point() +
  theme(text = element_text(size=15)) +
  facet_wrap(~maxHTO,scales="free_y") +
  ggtitle("Raw counts per maximum and second largest value HTO")


### DEMULTIPLEXING

HTO_tSNE_freemuxlet_type <- DimPlot(seuratObj, reduction = "HTO_tsne",
                                    group.by = "freemuxlet_type",
                                    pt.size = 1)
HTO_tSNE_freemuxlet_guess <- DimPlot(seuratObj,
                                     reduction = "HTO_tsne",
                                     group.by = "freemuxlet_guess",
                                     pt.size = 1)

secmax_DemuxLabel <- ggplot(HTO.data,aes(x = log2(maxval+1),y = log2(secval+1),color=Demux_label)) +
  geom_point(size=0.75) +
  geom_abline(slope = 1 ,color = "blue" ) +
  geom_abline(slope = 0.80 ,color = "red" ) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  facet_wrap(~maxHTO) +
  #scale_color_manual(values=c("firebrick", "coral","limegreen")) + #Doublet Negative Singlet
  ggtitle("Log2 counts per maximum and second largest value HTO: HTODemux")

secmax_MultiLabel <- ggplot(HTO.data,aes(x = log2(maxval+1),y = log2(secval+1),color=Multi_label)) +
  geom_point(size=0.75) +
  geom_abline(slope = 1 ,color = "blue" ) +
  geom_abline(slope = 0.80 ,color = "red" ) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  facet_wrap(~maxHTO) +
  #scale_color_manual(values=c("firebrick", "coral","limegreen")) + #Doublet Negative Singlet
  ggtitle("Log2 counts per maximum and second largest value HTO: Multiseq_Demux")

secmax_freemuxletLabel <- ggplot(HTO.data,aes(x = log2(maxval+1),y = log2(secval+1),color=freemuxlet_label)) +
  geom_point(size=0.75) +
  geom_abline(slope = 1 ,color = "blue" ) +
  geom_abline(slope = 0.80 ,color = "red" ) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  facet_wrap(~maxHTO) +
  #scale_color_manual(values=c("firebrick", "coral","limegreen")) + #Doublet Negative Singlet
  ggtitle("Log2 counts per maximum and second largest value HTO: Freemuxlet")

# HTO DEMUX
Demux_RNA_QQ <- ggplot(HTO.data, aes(x = nUMI, y = nGene, colour = Demux_class)) +
  geom_point(size = 2, aes(shape = Demux_class))+
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Unique Genes") +
  ggtitle(paste("Demux")) +
  theme_bw() +
  scale_color_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  theme(text = element_text(size=15))
Demux_RNA_nUMI <- ggplot(HTO.data, aes(x = Demux_class, y = nUMI)) +
  geom_violin(aes(fill=Demux_class)) +
  geom_jitter(size=0.2) +
  scale_fill_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  ylab("Count depth (total UMI count)") +
  ggtitle(paste("Demux")) +
  theme_bw() +
  theme(text = element_text(size=15))
Demux_RNA_nGene <- ggplot(HTO.data, aes(x = Demux_class, y = nGene)) +
  geom_violin(aes(fill=Demux_class)) +
  geom_jitter(size=0.2) +
  ylab("Nr of Unique Genes") +
  scale_fill_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  ggtitle(paste("Demux")) +
  theme_bw() +
  theme(text = element_text(size=15))
# MULTISEQ
Multi_RNA_QQ <- ggplot(HTO.data, aes(x = nUMI, y = nGene, colour = Multi_class)) +
  geom_point(size = 2, aes(shape = Multi_class))+
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Unique Genes") +
  scale_color_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  ggtitle(paste("Multi")) +
  theme_bw() +
  theme(text = element_text(size=15))
Multi_RNA_nUMI <- ggplot(HTO.data, aes(x = Multi_class, y = nUMI)) +
  geom_violin(aes(fill=Multi_class)) +
  geom_jitter(size=0.2) +
  ylab("Count depth (total UMI count)") +
  scale_fill_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  ggtitle(paste("Multi")) +
  theme_bw() +
  theme(text = element_text(size=15))
Multi_RNA_nGene <- ggplot(HTO.data, aes(x = Multi_class, y = nGene)) +
  geom_violin(aes(fill=Multi_class)) +
  geom_jitter(size=0.2) +
  scale_fill_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  ggtitle(paste("Multi")) +
  theme_bw() +
  theme(text = element_text(size=15))
# INHOUSE MIX
group.colors <- c(Singlet = "limegreen", Doublet = "firebrick", noNormalMixtures ="darkgrey")
mix_RNA_QQ <- ggplot(HTO.data, aes(x = nUMI, y = nGene, colour = mix.crit)) +
  geom_point(size = 2, aes(shape = mix.crit))+
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Unique Genes") +
  ggtitle(paste("Inhouse mixture")) +
  scale_color_manual(values=group.colors) + 
  theme_bw() +
  theme(text = element_text(size=15))
mix_RNA_nUMI <- ggplot(HTO.data, aes(x = mix.crit, y = nUMI)) +
  geom_violin(aes(fill=mix.crit)) +
  geom_jitter(size=0.2) +
  ylab("Count depth (total UMI count)") +
  ggtitle(paste("Inhouse mixture")) +
  theme_bw() +
  scale_fill_manual(values=group.colors) + 
  theme(text = element_text(size=15))
mix_RNA_nGene <- ggplot(HTO.data, aes(x = mix.crit, y = nGene)) +
  geom_violin(aes(fill=mix.crit)) +
  geom_jitter(size=0.2) +
  ggtitle(paste("Inhouse mixture")) +
  scale_fill_manual(values=group.colors) + 
  theme_bw() +
  theme(text = element_text(size=15))
# INHOUSE DIF
dif_RNA_QQ <- ggplot(HTO.data, aes(x = nUMI, y = nGene, colour = dif.crit)) +
  geom_point(size = 2, aes(shape = dif.crit))+
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Unique Genes") +
  ggtitle(paste("Inhouse difference")) +
  scale_color_manual(values=c("limegreen","firebrick")) + #Singlet Doublet
  theme_bw() +
  theme(text = element_text(size=15))
dif_RNA_nUMI <- ggplot(HTO.data, aes(x = dif.crit, y = nUMI)) +
  geom_violin(aes(fill=dif.crit)) +
  geom_jitter(size=0.2) +
  ylab("Count depth (total UMI count)") +
  ggtitle(paste("Inhouse difference")) +
  scale_fill_manual(values=c("limegreen","firebrick")) + #Singlet Doublet
  theme_bw() +
  theme(text = element_text(size=15))
dif_RNA_nGene <- ggplot(HTO.data, aes(x = dif.crit, y = nGene)) +
  geom_violin(aes(fill=dif.crit)) +
  geom_jitter(size=0.2) +
  ggtitle(paste("Inhouse difference")) +
  scale_fill_manual(values=c("limegreen","firebrick")) + #Singlet Doublet
  theme_bw() +
  theme(text = element_text(size=15))
# INHOUSE WITHIN
w_RNA_QQ <- ggplot(HTO.data, aes(x = nUMI, y = nGene, colour = w.crit)) +
  geom_point(size = 2, aes(shape = w.crit))+
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Unique Genes") +
  scale_color_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  ggtitle(paste("Inhouse within")) +
  theme_bw() +
  theme(text = element_text(size=15))
w_RNA_nUMI <- ggplot(HTO.data, aes(x = w.crit, y = nUMI)) +
  geom_violin(aes(fill=w.crit)) +
  geom_jitter(size=0.2) +
  scale_fill_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  ylab("Count depth (total UMI count)") +
  ggtitle(paste("Inhouse within")) +
  theme_bw() +
  theme(text = element_text(size=15))
w_RNA_nGene <- ggplot(HTO.data, aes(x = w.crit, y = nGene)) +
  geom_violin(aes(fill=w.crit)) +
  geom_jitter(size=0.2) +
  scale_fill_manual(values=c("limegreen","firebrick", "coral")) + #Singlet Doublet Negative
  ggtitle(paste("Inhouse within")) +
  theme_bw() +
  theme(text = element_text(size=15))
#COMBOS
# HTO_nUMI <- ggarrange(Demux_RNA_nUMI,Multi_RNA_nUMI,w_RNA_nUMI,
#                       dif_RNA_nUMI,ncol = 3,nrow=2)
#
# HTO_nGene <- ggarrange(Demux_RNA_nGene,Multi_RNA_nGene,w_RNA_nGene,
#                        dif_RNA_nGene,ncol = 3,nrow=2)
#
#
# HTO_QQ <- ggarrange(Demux_RNA_QQ,Multi_RNA_QQ,w_RNA_QQ,
#                     dif_RNA_QQ,ncol = 3,nrow=2)
HTO_nUMI <- ggarrange(Demux_RNA_nUMI,Multi_RNA_nUMI,w_RNA_nUMI,
                      mix_RNA_nUMI,dif_RNA_nUMI,ncol = 3,nrow=2)

HTO_nGene <- ggarrange(Demux_RNA_nGene,Multi_RNA_nGene,w_RNA_nGene,
                       mix_RNA_nGene,dif_RNA_nGene,ncol = 3,nrow=2)


HTO_QQ <- ggarrange(Demux_RNA_QQ,Multi_RNA_QQ,w_RNA_QQ,
                    mix_RNA_QQ,dif_RNA_QQ,ncol = 3,nrow=2)



############################
########## EXTRA PLOTS
############################
print("########################## extra visualisations ##########################")

counts <- as.matrix(GetAssayData(object = seuratObj,
                                 assay = "HTO",
                                 slot = "counts"))
ncells <- dim(counts)[2]
HashTags_present <- rownames(counts)
MULTI_ID <- HTO.data$Multi_label
MULTI_ID_df <- as.data.frame(MULTI_ID) %>%
  group_by(MULTI_ID) %>%  dplyr::count() %>% dplyr::mutate(perc = round(100*n/ncells,2))
HTODemux <- HTO.data$Demux_label
HTODemux_df <- as.data.frame(HTODemux) %>%
  group_by(HTODemux) %>%  dplyr::count() %>% dplyr::mutate(perc = round(100*n/ncells,2))

Inhouse_comp <- HTO.data$comp_label
Inhouse_comp_df <- as.data.frame(Inhouse_comp) %>%
  group_by(Inhouse_comp) %>%  dplyr::count() %>% dplyr::mutate(perc = round(100*n/ncells,2))
Inhouse_mix <- HTO.data$mix_label
Inhouse_mix_df <- as.data.frame(Inhouse_mix) %>%
  group_by(Inhouse_mix) %>%  dplyr::count() %>% dplyr::mutate(perc = round(100*n/ncells,2))

MAXHTO <- HTO.data$maxHTO
MAXHTO_df <- as.data.frame(MAXHTO) %>%
  group_by(MAXHTO) %>%  dplyr::count() %>% dplyr::mutate(perc = round(100*n/ncells,2))

MULTIseqDemux_results <- ggplot(MULTI_ID_df, aes(x = MULTI_ID, y = perc)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = paste("",perc,"%\nn=",n,"")), size=3.5, vjust = -0.3) +
  ggtitle("Classification according to MULTIseqDemux") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Percentage of cells") +
  ylim(0, 125) +
  theme(axis.text.x = element_text(angle = 45, size=12, hjust = 1))
HTODemux_results <- ggplot(HTODemux_df, aes(x = HTODemux, y = perc)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = paste("",perc,"%\nn=",n,"")), size=3.5, vjust = -0.3) +
  ggtitle("Classification according to HTODemux") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Percentage of cells") +
  ylim(0, 125) +
  theme(axis.text.x = element_text(angle = 45, size=12, hjust = 1))

Inhouse_comp_results <- ggplot(Inhouse_comp_df, aes(x = Inhouse_comp, y = perc)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = paste("",perc,"%\nn=",n,"")), size=3.5, vjust = -0.3) +
  ggtitle("Classification according to Inhouse_comp") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Percentage of cells") +
  ylim(0, 125) +
  theme(axis.text.x = element_text(angle = 45, size=12, hjust = 1))
Inhouse_mix_results <- ggplot(Inhouse_mix_df, aes(x = Inhouse_mix, y = perc)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = paste("",perc,"%\nn=",n,"")), size=3.5, vjust = -0.3) +
  ggtitle("Classification according to Inhouse_mix") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Percentage of cells") +
  ylim(0, 125) +
  theme(axis.text.x = element_text(angle = 45, size=12, hjust = 1))

MAXHTO_results <- ggplot(MAXHTO_df, aes(x = MAXHTO, y = perc)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = paste("",perc,"%\nn=",n,"")), size=3.5, vjust = -0.3) +
  ggtitle("Cell frequency per maximum hashtag") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Percentage of cells") +
  ylim(0, 125) +
  theme(axis.text.x = element_text(angle = 45, size=12, hjust = 1))

FREEMUXLET_ID <- seuratObj$freemuxlet_label

FREEMUXLET_ID_df <- as.data.frame(FREEMUXLET_ID) %>%
  group_by(FREEMUXLET_ID) %>%  
  dplyr::count() %>% 
  dplyr::mutate(perc = round(100*n/ncells,2))

FREEMUXLET_DEMUX_results <- ggplot(FREEMUXLET_ID_df, aes(x = FREEMUXLET_ID, y = perc)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = paste("",perc,"%\nn=",n,"")), size=3.5, vjust = -0.3) +
  ggtitle("Classification according to MULTIseqDemux") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Percentage of cells") +
  ylim(0, 125) +
  theme(axis.text.x = element_text(angle = 45, size=12, hjust = 1))


############################
########## SAVE HTO PLOTS
############################

print("########################## saving plots ##########################")
# png(filename = paste0(output.dir,"Plots/HTO/","01_HashDistribution"), width = 1100, height = 700)
# max.barplot
# dev.off()

#HTO_QC1<- (MAXHTO_results + HTODemux_results)/MULTIseqDemux_results

# HASH DISTRIBUTION
png(filename = paste0("Plots/HTO/","01_HashDistribution_freemuxlet",".png"), width = 1000, height = 1000, type = 'cairo')
FREEMUXLET_DEMUX_results 
dev.off()

png(filename = paste0("Plots/HTO/","01_HashDistribution_maxHTO",".png"), width = 700, height = 700, type = 'cairo')
MAXHTO_results
dev.off()

png(filename = paste0("Plots/HTO/","01_HashDistribution_demultiplexing",".png"),
    width = 1400,
    height = 700,
    type = 'cairo')
ggarrange(HTODemux_results,MULTIseqDemux_results,ncol=2)
dev.off()

png(filename = paste0("Plots/HTO/","01_HashDistribution_inhouse",".png"), width = 700, height = 700, type = 'cairo')
Inhouse_comp_results
dev.off()

png(filename = paste0("Plots/HTO/","01_HashDistribution",".png"), width = 1500, height = 1000, type = 'cairo')
ggarrange(HTODemux_results,MULTIseqDemux_results,
          MAXHTO_results,
          Inhouse_comp_results,
          ncol=2, nrow=2)
dev.off()

pdf(file = paste0("Plots/HTO/","01_HashDistribution",".pdf"), width = 10, height = 10)
HTODemux_results
MULTIseqDemux_results
# Inhouse_comp_results
MAXHTO_results
ggarrange(HTODemux_results,MULTIseqDemux_results,
          MAXHTO_results,
          Inhouse_comp_results,
          ncol=2, nrow=2)
dev.off()

# HASH DENSITY
png(filename = paste0("Plots/HTO/","02_HashDensity",".png"), width = 1100, height = 700, type = 'cairo')
ggarrange(HTO_dens, secmax_dens, ncol = 1)
dev.off()
pdf(file = paste0("Plots/HTO/","02_HashDensity",".pdf"), width = 11, height = 7)
ggarrange(HTO_dens, secmax_dens, ncol = 1)
dev.off()
# DOUBLET CHECK
png(filename = paste0("Plots/HTO/","03b_Doubletcheck_HTODemux",".png"), width = 500, height = 500, type = 'cairo')
secmax_DemuxLabel
dev.off()
png(filename = paste0("Plots/HTO/","03a_Doubletcheck_Inhouse",".png"), width = 500, height = 500, type = 'cairo')
secmax_mix
dev.off()
png(filename = paste0("Plots/HTO/","03c_Doubletcheck_MultiSeq",".png"), width = 500, height = 500, type = 'cairo')
secmax_MultiLabel
dev.off()
png(filename = paste0("Plots/HTO/","03c_Doubletcheck_Freemuxlet",".png"), width = 500, height = 500, type = 'cairo')
secmax_freemuxletLabel
dev.off()
png(filename = paste0("Plots/HTO/","03_Doubletcheck",".png"), width = 1000, height = 1000, type = 'cairo')
ggarrange(secmax_mix,secmax_DemuxLabel,secmax_MultiLabel,secmax_freemuxletLabel, nrow = 2, ncol = 2)
dev.off()
pdf(file = paste0("Plots/HTO/","03_Doubletcheck",".pdf"), width = 20, height = 10)
ggarrange(secmax_mix,secmax_DemuxLabel,secmax_MultiLabel,secmax_freemuxletLabel, nrow = 2, ncol = 2)
dev.off()
# COUNT DISTRIBUTION
png(filename = paste0("Plots/HTO/","04_CountDistribution",".png"), width = 1100, height = 700, type = 'cairo')
HTO_cnt
dev.off()
pdf(file = paste0("Plots/HTO/","04_CountDistribution",".pdf"), width = 11, height = 7)
HTO_cnt
dev.off()
# RAW COUNTS
png(filename = paste0("Plots/HTO/","05a_RawCounts",".png"), width = 1100, height = 700, type = 'cairo')
HTO_dif2
dev.off()
png(filename = paste0("Plots/HTO/","05b_RawCounts",".png"), width = 1100, height = 700, type = 'cairo')
HTO_lines
dev.off()
pdf(file = "Plots/HTO/05_RawCounts.pdf", width = 11, height = 7)
HTO_dif2
HTO_lines
dev.off()
# VIOLIN
png(filename = "Plots/HTO/06_Violin_UMI.png", width = 1500, height = 1000, type = 'cairo')
HTO_nUMI
dev.off()
png(filename = "Plots/HTO/06_Violin_nGene.png", width = 1500, height = 1000, type = 'cairo')
HTO_nGene
dev.off()
png(filename = "Plots/HTO/06_Violin_QQ.png", width = 1500, height = 1000, type = 'cairo')
HTO_QQ
dev.off()
pdf(file = "Plots/HTO/06_Violin.pdf", width=15, height = 10)
HTO_nUMI
HTO_nGene
HTO_QQ
dev.off()
# TSNE
png(filename = "Plots/HTO/07_tSNE_HTO_HTODemux.ID.png", width = 1000, height = 1000, type = 'cairo')
DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "hash.ID", pt.size = 1)
dev.off()
png(filename = "Plots/HTO/07_tSNE_HTO_HTODemux.png", width = 1000, height = 1000, type = 'cairo')
DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "HTO_classification.global", pt.size = 1)
dev.off()
png(filename = "Plots/HTO/07_tSNE_HTO_multiID.png", width = 1000, height = 1000, type = 'cairo')
DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "MULTI_ID", pt.size = 1)
dev.off()
png(filename = "Plots/HTO/07_tSNE_HTO_compare.png", width = 1500, height = 500, type = 'cairo')
ggarrange(DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "hash.ID", pt.size = 1),
          DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "HTO_classification.global", pt.size = 1),
          DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "MULTI_ID", pt.size = 1),
          ncol=3)
dev.off()
pdf(file = "Plots/HTO/07_tSNE_HTO.pdf", width=15, height = 10)
DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "hash.ID", pt.size = 1)
DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "HTO_classification.global", pt.size = 1)
DimPlot(seuratObj, reduction = "HTO_tsne", group.by = "MULTI_ID", pt.size = 1)
dev.off()
png(filename = "Plots/HTO/08_tSNE_HTO_freemuxlet.png", width = 1400, height = 700, type = 'cairo')
ggarrange(HTO_tSNE_freemuxlet_type,HTO_tSNE_freemuxlet_guess,ncol=2)
dev.off()
png(filename = "Plots/HTO/08_tSNE_HTO_compare.png", width = 1500, height = 1000, type = 'cairo')
ggarrange(DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "freemuxlet_guess",
                  pt.size = 1) +
            ggtitle("Freemuxlet:label"),
          DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "hash.ID",
                  pt.size = 1) +
            ggtitle("HTOdemux:label"),
          DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "MULTI_ID",
                  pt.size = 1) +
            ggtitle("Multiseq:label"),
          DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "freemuxlet_type",
                  pt.size = 1) +
            ggtitle("Freemuxlet:type"),
          DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "HTO_classification.global",
                  pt.size = 1) +
            ggtitle("HTOdemux:type"),
          ncol=3, nrow=2)
dev.off()
pdf(file =  "Plots/HTO/08_tSNE_HTO.pdf", width=10, height = 10)
DimPlot(seuratObj,
        reduction = "HTO_tsne",
        group.by = "freemuxlet_guess",
        pt.size = 1) +
  ggtitle("Freemuxlet:label")
DimPlot(seuratObj,
        reduction = "HTO_tsne",
        group.by = "hash.ID",
        pt.size = 1) +
  ggtitle("HTOdemux:label")
DimPlot(seuratObj,
        reduction = "HTO_tsne",
        group.by = "MULTI_ID",
        pt.size = 1) +
  ggtitle("Multiseq:label")
DimPlot(seuratObj,
        reduction = "HTO_tsne",
        group.by = "freemuxlet_type",
        pt.size = 1) +
  ggtitle("Freemuxlet:type")
DimPlot(seuratObj,
        reduction = "HTO_tsne",
        group.by = "HTO_classification.global",
        pt.size = 1) +
  ggtitle("HTOdemux:type")
dev.off()

png(filename = "Plots/HTO/09_tSNE_HTO_compare.png",
    width = 1500,
    height = 1000,
    type = 'cairo')
ggarrange(DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "freemuxlet_guess",
                  pt.size = 1,
                  cols=length(unique(seuratObj$freemuxlet_guess))) +
            ggtitle("Freemuxlet:label"),
          DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "hash.ID",
                  pt.size = 1,
                  cols=length(unique(seuratObj$hash.ID))) +
            ggtitle("HTOdemux:label"),
          DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "MULTI_ID",
                  pt.size = 1,
                  cols=length(unique(seuratObj$MULTI_ID))) +
            ggtitle("Multiseq:label"),
          DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "freemuxlet_type",
                  pt.size = 1,
                  cols=length(unique(seuratObj$freemuxlet_type))) +
            ggtitle("Freemuxlet:type"),
          DimPlot(seuratObj,
                  reduction = "HTO_tsne",
                  group.by = "HTO_classification.global",
                  pt.size = 1,
                  cols=length(unique(seuratObj$HTO_classification.global))) +
            ggtitle("HTOdemux:type"),
          ncol=3, nrow=2)
dev.off()
png(filename = "Plots/HTO/09_confusionmatrix.png",
    width = 500,
    height = 300,
    type = 'cairo')
ggplot(data = TheMatrix, mapping = aes(x = Reference, y = Prediction, fill="White")) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(Singlet = "green", bad = "red")) +
  theme_bw() +
  xlim(rev(levels(TheMatrix$Reference))) +
  xlab("HTOdemux_globalclass") + ylab(label = "Freemuxlet")
dev.off()

pdf(file = "Plots/HTO/Overview_COV010.pdf", width = 10, height = 10)

lapply(seuratObj@misc$PresentHashingAntibodies,function(x) {
  FeaturePlot(seuratObj,features = x, pt.size=0.10,label.size = 8,order=T)
})

FeaturePlot(seuratObj,
            features=(seuratObj@misc$PresentHashingAntibodies),
            ncol=2,
            pt.size=0.1,
            label.size = 8,
            order=T)

DimPlot(seuratObj,
        group.by = "MULTI_ID",
        label = F,
        pt.size = 0.1) +
  ggtitle("MULTIseqDemux")
DimPlot(seuratObj,
        split.by = "MULTI_ID",
        group.by = "MULTI_ID",
        label = F,
        pt.size = 0.1,
        ncol = 2) +
  ggtitle("MULTIseqDemux")
DimPlot(seuratObj,
        group.by = "HTO_maxID",
        label = F,
        pt.size = 0.1) +
  ggtitle("MAX HTO count")
DimPlot(seuratObj,
        split.by = "HTO_maxID",
        group.by = "HTO_maxID",
        label = F,
        pt.size = 0.1,
        ncol = 2) +
  ggtitle("MAX HTO count")
DimPlot(seuratObj,
        group.by = "freemuxlet_guess",
        label = F,
        pt.size = 0.1) +
  ggtitle("freemuxlet")
DimPlot(seuratObj,
        split.by = "freemuxlet_guess",
        group.by = "freemuxlet_guess",
        label = F,
        pt.size = 0.1,
        ncol = 2) +
  ggtitle("freemuxlet")
DimPlot(seuratObj,
        group.by = "subsample",
        label = F,
        pt.size = 0.1) +
  ggtitle("subsamples")
DimPlot(seuratObj,
        split.by = "subsample",
        group.by = "subsample",
        label = F,
        pt.size = 0.1,
        ncol = 2) +
  ggtitle("subsamples")
dev.off()

lapply(seuratObj@misc$subsamples,function(x){
  if(sum(seuratObj$subsample==x) > 0 ){
    saveRDS(subset(seuratObj, subset = subsample == x),file = paste0("Robjects/seuratObj_",seuratObj@project.name,"_hash",strsplit(x,".#")[[1]][2],".rds"))
  }else{
    print(paste0("No cell remaining for the subsample ",x))
  }
})

saveRDS(seuratObj, file = paste0(seuratObj@project.name,".SEURAT__HTO_VISUALISATION.rds"))