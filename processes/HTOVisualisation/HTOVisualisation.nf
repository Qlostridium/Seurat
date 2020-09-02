#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process SEURAT__HTO_VISUALISATION {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	publishDir "${params.folder.outDir}/${samplename}", mode: 'symlink', pattern : "Robjects/**"
	publishDir "${params.folder.outDir}/${samplename}", mode: 'move', pattern: "Plots/HTO/**"
	container params.Seurat.container
  input:
	tuple val(samplename), file(seuratobj)
  output:
	tuple val(samplename), file("${samplename}.SEURAT__HTO_VISUALISATION.rds")
	file("Robjects/**")
	file("Plots/HTO/**")
  script:
	"""
	Rscript ${scriptDir}/HTOVisualisation/HTOVisualisation.R --seuratObj ${seuratobj} \
	--scriptFunctions ${scriptDir}/utils/script_functions_COVID.R
	"""
}
