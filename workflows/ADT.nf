#!/usr/bin/env nextflow
nextflow.preview.dsl=2
include SEURAT__NORMALIZATION from '../processes/normalization/normalization.nf' params(params)
include SEURAT__SCALING from '../processes/scaling/scaling.nf' params(params)
include SEURAT__FIND_VARIABLE_FEATURES from '../processes/findVariableFeatures/findVariableFeatures.nf' params(params)
include SEURAT__DIMENSIONALITY_REDUCTION_PCA from '../processes/dimensionalityReduction/dimensionalityReduction.nf' params(params)
include SEURAT__CLUSTERING from '../processes/clustering/clustering.nf' params(params)

workflow run_ADT {
	take: inputtuple
	main :
		SEURAT__NORMALIZATION(inputtuple,"ADT")
		SEURAT__SCALING(SEURAT__NORMALIZATION.out,"ADT")
		SEURAT__FIND_VARIABLE_FEATURES(SEURAT__SCALING.out,"ADT")
		SEURAT__DIMENSIONALITY_REDUCTION_PCA(SEURAT__FIND_VARIABLE_FEATURES.out[0],"ADT")
		SEURAT__CLUSTERING(SEURAT__DIMENSIONALITY_REDUCTION_PCA.out[0],"ADT")
}
