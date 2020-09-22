#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__THRESHOLDFILTERING{
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink', pattern: "Plots/RNA/**"
	input:
	tuple val(samplename), file(sobj)
	output:
	tuple val(samplename), file("${samplename}.SEURAT__THRESHOLDFILTERING.rds")
	file("Plots/RNA/**")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.thresholdFiltering)
	"""
	Rscript ${scriptDir}/thresholdFiltering/thresholdFiltering.R --seuratObj "${sobj}" \
		--output "${samplename}.SEURAT__THRESHOLDFILTERING.rds" \
		--scriptFunctions ${scriptDir}/utils/script_functions_COVID.R \
		${(!(sampleParams.nmad_low_feature == null)) ? '--nmad_low_feature ' + sampleParams.nmad_low_feature: ''} \
		${(sampleParams.nmad_high_feature == null) ? '' : '--nmad_high_feature ' + sampleParams.nmad_high_feature} \
		${(sampleParams.nmad_low_UMI == null) ? '' : '--nmad_low_UMI ' + sampleParams.nmad_low_UMI} \
		${(sampleParams.nmad_high_UMI == null) ? '' : '--nmad_high_UMI ' + sampleParams.nmad_high_UMI} \
		${(sampleParams.nmad_high_mito == null) ? '' : '--nmad_high_mito ' + sampleParams.nmad_high_mito} \

	"""
}
