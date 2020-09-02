#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process SEURAT__SCE_TO_SEURAT_WITH_MERGE {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	container params.Seurat.container
  input:
	tuple val(samplename), file(sceobj)
	tuple val(samplename), file(seuratobj)
  output:
	tuple val(samplename),file("${samplename}.SEURAT__SCE_TO_SEURAT_WITH_MERGE.rds")
  script:
	"""
	Rscript ${scriptDir}/utils/sceToSeurat.R --sceObj ${sceobj} \
	--seuratObj ${seuratobj} \
	--output "${samplename}.SEURAT__SCE_TO_SEURAT_WITH_MERGE.rds"
	"""
}

process SEURAT__SCE_TO_SEURAT {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	container params.Seurat.container
  input:
	tuple val(samplename), file(sceobj)
  output:
	tuple val(samplename),file("${samplename}.SEURAT__SCE_TO_SEURAT.rds") \
	--output "${samplename}.SEURAT__SCE_TO_SEURAT.rds"
  script:
	"""
	Rscript ${scriptDir}/utils/sceToSeurat.R --sceObj ${sceobj}
	"""
}

process SEURAT__SEURAT_TO_SCE {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	container params.Seurat.container
  input:
	tuple val(samplename), file(seuratobj)
  output:
	tuple val(samplename),file("${samplename}.SEURAT__SEURAT_TO_SCE.rds")
  script:
	"""
	Rscript ${scriptDir}/utils/seuratToSce.R --seuratObj ${seuratobj} \
	--output "${samplename}.SEURAT__SEURAT_TO_SCE.rds"
	"""
}