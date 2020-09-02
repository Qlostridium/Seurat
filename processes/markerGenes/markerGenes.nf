#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process SEURAT__MARKER_GENES{
	publishDir "${params.folder.outDir}/${samplename}", mode: 'symlink'
	input:
	tuple val(samplename), file(sobj)
	output:
	tuple val(samplename), file("${samplename}.SEURAT__MARKER_GENES.rds")
	file("Plots/**")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.markerGenes)
	"""
	Rscript ${scriptDir}/markerGenes/markerGenes.R --seuratObj "${sobj}" \
		--output "${samplename}.SEURAT__MARKER_GENES.rds" \
		--markerGenes "${sampleParams.markerGenes}"
	"""
}
