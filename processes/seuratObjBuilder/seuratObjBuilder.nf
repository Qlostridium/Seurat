#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process SEURAT__SEURAT_OBJECT_BUILDER {
	//publishDir "${params.folder.outDir}/${samplename}", mode: 'symlink'
	container params.Seurat.container
	input:
	tuple val(samplename), val(cmat)
	output:
	tuple val(samplename), file("${samplename}.SEURAT__SEURAT_OBJECT_BUILDER.rds")
	script:
	"""
	Rscript ${scriptDir}/seuratObjBuilder/seuratObjBuilder.R --inputMatrix "${cmat}" \
		--output "${samplename}.SEURAT__SEURAT_OBJECT_BUILDER.rds" \
		--sample $samplename \
		${(!(params.Seurat.seuratObjBuilder.minFeaturesGEX == null)) ? '--minFeaturesGEX ' + params.Seurat.seuratObjBuilder.minFeaturesGEX: ''} \
		${(!(params.Seurat.seuratObjBuilder.minCellsGEX == null)) ? '--minCellsGEX ' + params.Seurat.seuratObjBuilder.minCellsGEX: ''}
	"""
}
