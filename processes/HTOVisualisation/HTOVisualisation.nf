#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (!(params.folder.standAlone == null) && params.folder.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"
process SEURAT__HTO_VISUALISATION {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	publishDir "${params.global.outDir}/${samplename}", mode: 'symlink', pattern : "Robjects/**"
	publishDir "${params.global.outDir}/${samplename}", mode: 'move', pattern: "Plots/HTO/**"
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
