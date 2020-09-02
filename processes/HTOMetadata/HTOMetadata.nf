#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process SEURAT__HTO_METADATA {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	container params.Seurat.container
  input:
	tuple val(samplename), file(seuratobj)
  output:
	tuple val(samplename), file("${samplename}.SEURAT__HTO_METADATA.rds")
  script:
	"""
	Rscript ${scriptDir}/HTOMetadata/HTOMetadata.R --seuratObj ${seuratobj} \
	--output "${samplename}.SEURAT__HTO_METADATA.rds" \
	--scriptFunctions ${scriptDir}/utils/script_functions_COVID.R
	"""
}
