#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process SEURAT__RNA_QC{
	publishDir "${params.folder.outDir}/${samplename}", mode: 'move', pattern: "Plots/RNA/**"
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
