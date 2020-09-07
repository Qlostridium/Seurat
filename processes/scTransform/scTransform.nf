#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (!(params.folder.standAlone == null) && params.folder.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__SCTRANSFORM{
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink'
	input:
	tuple val(samplename), file(sobj)
	output:
	tuple val(samplename), file("${samplename}.SEURAT__SCTRANSFORM.rds")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.scTransform)
	"""
	Rscript ${scriptDir}/scTransform/scTransform.R --seuratObj "${sobj}" \
		--output "${samplename}.SEURAT__SCTRANSFORM.rds" \
		${(sampleParams.regressSubsamples == null || sampleParams.regressSubsamples == "false") ? '': '--regressSubsamples'}
	"""
}
