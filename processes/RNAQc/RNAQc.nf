#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__RNA_QC{
	publishDir "${params.global.outdir}/${samplename}", mode: 'move', pattern: "Plots/RNA/**"
	container params.Seurat.container
	input:
	tuple val(samplename), file(sobj)
	output:
	tuple val(samplename), file("${samplename}.SEURAT__RNA_QC.rds")
	file("Plots/RNA/**")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.RNAQc)
	"""
	Rscript ${scriptDir}/RNAQc/RNAQc.R --seuratObj "${sobj}" \
		--output "${samplename}.SEURAT__RNA_QC.rds" \
		${(sampleParams.mitoGenes == null || sampleParams.mitoGenes == "false") ? '' : '--mitoGenes'} \
		${(sampleParams.covidGenes == null || sampleParams.covidGenes == "false") ? '' : '--covidGenes'} \
		${(sampleParams.rbcGenes == null || sampleParams.rbcGenes == "false") ? '' : '--rbcGenes'} \
		${(sampleParams.genomeName1 == null) ? '' : '--genomeName1 ' + sampleParams.genomeName1} \
		${(sampleParams.genomeName2 == null) ? '' : '--genomeName2 ' + sampleParams.genomeName2}
	"""
}
