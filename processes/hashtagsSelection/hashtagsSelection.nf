#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process SEURAT__HASHTAGS_SELECTION{
	publishDir "${params.folder.outDir}/${samplename}", mode: 'symlink', pattern : "${samplename}_logQC.txt"
	container params.Seurat.container
	input:
	tuple val(samplename), file(sobj)
	output:
	tuple val(samplename), file("${samplename}.SEURAT__HASHTAG_SELECTION.rds")
	file("${samplename}_logQC.txt")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.hashtagsSelection)
	"""
	Rscript ${scriptDir}/hashtagsSelection/hashtagsSelection.R --inputSeuratRds "${sobj}" \
		--output "${samplename}.SEURAT__HASHTAG_SELECTION.rds" \
		--hashtags "${sampleParams.hashtags}"
	"""
}
