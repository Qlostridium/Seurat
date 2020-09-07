#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (!(params.global.standAlone == null) && params.global.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

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
