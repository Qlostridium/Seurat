#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (!(params.global.standAlone == null) && params.global.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__NORMALIZATION{
	//publishDir "${params.global.outDir}/${samplename}", mode: 'symlink'
	container params.Seurat.container
	input:
	tuple val(samplename), file(sobj)
	val assayType
	output:
	tuple val(samplename), file("${samplename}.SEURAT__NORMALIZATION_${assayType}.rds")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.normalization)
	def assayParams = params.configParser(assayType, sampleParams)
	"""
	Rscript ${scriptDir}/normalization/normalization.R --inputSeuratRds "${sobj}" \
		--output "${samplename}.SEURAT__NORMALIZATION_${assayType}.rds" \
		--assay "${assayType}" \
		${(!(assayParams.containsKey == null)) ? '--margin ' + assayParams.margin: ''} \
		${(!(assayParams.scalefactor == null)) ? '--scalefactor ' + assayParams.scalefactor: ''} \
		${(!(assayParams.normalizationMethod == null)) ? '--normalizationMethod ' + assayParams.normalizationMethod: ''}
	"""
}
