#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__FIND_VARIABLE_FEATURES{
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink', pattern: "Plots/**"
	container params.Seurat.container
	input:
	tuple val(samplename), file(sobj)
	val assayType
	output:
	tuple val(samplename), file("${samplename}.SEURAT__FIND_VARIABLE_FEATURES_${assayType}.rds")
	file("Plots/**")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.findVariableFeatures)
	def assayParams = params.configParser(assayType, sampleParams)
	"""
	Rscript ${scriptDir}/findVariableFeatures/findVariableFeatures.R --seuratObj "${sobj}" \
		--output "${samplename}.SEURAT__FIND_VARIABLE_FEATURES_${assayType}.rds" \
		--assay "$assayType" \
		${(!(assayParams.selectionMethod == null)) ? '--selectionMethod ' + assayParams.selectionMethod: ''} \
		${(!(assayParams.loesSpan == null)) ? '--loesSpan ' + assayParams.loesSpan: ''} \
		${(!(assayParams.clipMax == null)) ? '--clipMax ' + assayParams.clipMax: ''} \
		${(!(assayParams.meanFunction == null)) ? '--meanFunction ' + assayParams.meanFunction: ''} \
		${(!(assayParams.dispersionFunction == null)) ? '--dispersionFunction ' + assayParams.dispersionFunction: ''} \
		${(!(assayParams.numBin == null)) ? '--numBin ' + assayParams.numBin: ''} \
		${(!(assayParams.binningMethod == null)) ? '--binningMethod ' + assayParams.binningMethod: ''} \
		${(!(assayParams.nfeatures == null)) ? '--nfeatures ' + assayParams.nfeatures: ''} \
		${(!(assayParams.meanCutoff == null)) ? '--meanCutoff ' + assayParams.meanCutoff: ''} \
		${(!(assayParams.dispersionCutoff == null)) ? '--dispersionCutoff ' + assayParams.dispersionCutoff: ''} \
		${(!(assayParams.verbose == null)) ? '--verbose ' + assayParams.verbose: ''}
	"""
}
