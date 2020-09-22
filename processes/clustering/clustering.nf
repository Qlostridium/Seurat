#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__CLUSTERING {
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink', pattern: "Plots/**"
	input:
	tuple val(samplename), file(sobj)
	val assay
	output:
	tuple val(samplename), file("${samplename}.SEURAT__CLUSTERING_${assay}.rds")
	file("Plots/**")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.clustering)
	def assayParams = params.configParser(assay, sampleParams)
	"""
	Rscript ${scriptDir}/clustering/clustering.R --seuratObj "${sobj}" \
		--output "${samplename}.SEURAT__CLUSTERING_${assay}.rds" \
		--assay ${assay} \
		${(assayParams.dimsToUse == null) ? '' :'--dimsToUse ' +assayParams.dimsToUse } \
		${(assayParams.resToUse == null) ? '' :'--resToUse ' + assayParams.resToUse} \
		${(assayParams.perplexity == null) ? '' : '--perplexity ' + assayParams.perplexity} \
		${(assayParams.assayForCrossModalityGraphs == null) ? '' : '--assayForCrossModalityGraphs ' + assayParams.assayForCrossModalityGraphs}
	"""
}
