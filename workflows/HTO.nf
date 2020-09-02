#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include SEURAT__SEURAT_OBJECT_BUILDER from '../processes/seuratObjBuilder/seuratObjBuilder.nf' params(params)
include SEURAT__NORMALIZATION from '../processes/normalization/normalization.nf' params(params)
include SEURAT__HASHTAGS_SELECTION from '../processes/hashtagsSelection/hashtagsSelection.nf' params(params)
include SEURAT__DIMENSIONALITY_REDUCTION_TSNE from '../processes/dimensionalityReduction/dimensionalityReduction.nf' params(params)
include SEURAT__HTO_METADATA from '../processes/HTOMetadata/HTOMetadata.nf' params(params)
include SEURAT__HTO_DEMULTIPLEXNG from '../processes/HTODemultiplexing/HTODemultiplexing.nf' params(params)
include SEURAT__HTO_VISUALISATION from '../processes/HTOVisualisation/HTOVisualisation.nf' params(params)

workflow run_HTO {
	take: data
	main:
		SEURAT__SEURAT_OBJECT_BUILDER(data)
		SEURAT__HASHTAGS_SELECTION(SEURAT__SEURAT_OBJECT_BUILDER.out)
		SEURAT__NORMALIZATION(SEURAT__HASHTAGS_SELECTION.out[0],"HTO")
		SEURAT__DIMENSIONALITY_REDUCTION_TSNE(SEURAT__NORMALIZATION.out,"HTO")
		SEURAT__HTO_METADATA(SEURAT__DIMENSIONALITY_REDUCTION_TSNE.out)
		SEURAT__HTO_DEMULTIPLEXNG(SEURAT__HTO_METADATA.out,file(params.Seurat.HTODemultiplexing.freemuxletFile),SEURAT__HASHTAGS_SELECTION.out[1])
		SEURAT__HTO_VISUALISATION(SEURAT__HTO_DEMULTIPLEXNG.out[0])
	emit:
		SEURAT__HTO_VISUALISATION.out[0]
}
