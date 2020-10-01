#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__DIMENSIONALITY_REDUCTION_PCA {
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink', pattern: "Plots/**"
	container params.Seurat.container
	input:
	tuple val(samplename), file(seuratobj)
	val assay
	output:
	tuple val(samplename), file("${samplename}.SEURAT__DIMENSIONALITY_REDUCTION_PCA_${assay}.rds")
	file("Plots/**")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.dimensionalityReduction.pca)
	def assayParams = params.configParser(assay, sampleParams)
	"""
	Rscript ${scriptDir}/dimensionalityReduction/pca.R --seuratObj ${seuratobj} \
		--output "${samplename}.SEURAT__DIMENSIONALITY_REDUCTION_PCA_${assay}.rds" \
		--assay $assay \
		${(assayParams.nPcs == null) ? '' : '--nPcs ' + assayParams.nPcs} \
		${(assayParams.nPlotedPcs == null) ? '' : '--nPlotedPcs ' + assayParams.nPlotedPcs} \
		${(assayParams.scaleData == null || assayParams.scaleData != "true") ? '' :'--scaleData'} \
		${(assayParams.removeScaledData == null || assayParams.removeScaledData != "true") ? '' :'--removeScaledData'} \
		${(assayParams.excludePatternHVG == null) ? '' : '--excludePatternHVG ' + assayParams.excludePatternHVG} \
		${(assayParams.diagnosticPlots == null) ? '' : '--diagnosticPlots'}
	"""
}

process SEURAT__DIMENSIONALITY_REDUCTION_TSNE {
	//publishDir "${params.out_dir}/${samplename}", mode: 'copy'
	container params.Seurat.container
	input:
	tuple val(samplename), file(seuratobj)
	val assayType
	output:
	tuple val(samplename), file("${samplename}.SEURAT__DIMENSIONALITY_REDUCTION_TSNE_${assayType}.rds")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.dimensionalityReduction.tsne)
	def assayParams = params.configParser(assayType, sampleParams)
	"""
	Rscript ${scriptDir}/dimensionalityReduction/tsne.R --inputSeuratRds ${seuratobj} \
		--output "${samplename}.SEURAT__DIMENSIONALITY_REDUCTION_TSNE_${assayType}.rds" \
		--assay "${assayType}" \
		${(!(assayParams.dims == null)) ? '--dims ' + assayParams.dims: ''} \
		${(!(assayParams.reduction == null)) ? '--reduction ' + assayParams.reduction: ''} \
		${(!(assayParams.seedUse == null)) ? '--seedUse ' + assayParams.seedUse: ''} \
		${(!(assayParams.tsneMethod == null)) ? '--tsneMethod ' + assayParams.tsneMethod: ''} \
		${(!(assayParams.addIter == null)) ? '--addIter ' + assayParams.addIter: ''} \
		${(!(assayParams.dimEmbed == null)) ? '--dimEmbed ' + assayParams.dimEmbed: ''} \
		${(!(assayParams.perplexity == null)) ? '--perplexity ' + assayParams.perplexity: ''}
	"""
}

/*
process SEURAT__DIMANSIONALITY_REDUCTION_UMAP {
	//publishDir "${params.out_dir}/${samplename}", mode: 'copy'
	input:
	tuple val(samplename), file(seuratobj)
	output:
	tuple val(samplename), file("${samplename}.SEURAT__DIMANSIONALITY_REDUCTION_UMAP.rds")
	script:
	"""
	Rscript ${params.pdir}/modules/dimensionalityReduction/umap.R --inputSeuratRds ${seuratobj} \
		--output "${samplename}.SEURAT__DIMANSIONALITY_REDUCTION_UMAP.rds"
	"""
}
*/
