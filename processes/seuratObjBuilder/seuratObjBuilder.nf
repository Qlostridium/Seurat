#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (!(params.folder.standAlone == null) && params.folder.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__SEURAT_OBJECT_BUILDER {
	//publishDir "${params.folder.outDir}/${samplename}", mode: 'symlink'
	container params.Seurat.container
	input:
	tuple val(samplename), val(cmat)
	output:
	tuple val(samplename), file("${samplename}.SEURAT__SEURAT_OBJECT_BUILDER.rds")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.seuratObjBuilder)
	"""
	Rscript ${scriptDir}/seuratObjBuilder/seuratObjBuilder.R --inputMatrix "${cmat}" \
		--output "${samplename}.SEURAT__SEURAT_OBJECT_BUILDER.rds" \
		--sample $samplename \
		${(!(sampleParams.minFeaturesGEX == null)) ? '--minFeaturesGEX ' + sampleParams.minFeaturesGEX: ''} \
		${(!(sampleParams.minCellsGEX == null)) ? '--minCellsGEX ' + sampleParams.minCellsGEX: ''}
	"""
}
