#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone != true) ? "${params.global.rundir}/src/Seurat/processes" : "${params.global.rundir}/processes"

process annotation_graphs {
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink', pattern: "Plots/**"
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink', pattern: "allClusters_${samplename}.xlsx"
	container params.Seurat.container
  input:
	tuple val(samplename), file(seuratobj)
	file(markersfile)
	val(assay)
  output:
	file("Plots/***")
	file("allClusters_${samplename}.xlsx")
  script:
  	def realmarkersfile = markersfile.name != 'NO_FILE' ? "--markersFile $markersfile" : ''
	"""
	Rscript ${scriptDir}/annotationGraphs/annotationGraphs.R --seuratObj ${seuratobj} \
	--assay $assay \
	$realmarkersfile
	"""
}

workflow SEURAT__ANNOTATION_GRAPHS {
	take:
		input
		assayname
	main:

		if(params.Seurat.annotationGraphs.containsKey(assayname)){
			assayParams = params.Seurat.annotationGraphs."${assayname}"
		} else {
			assayParams = params.Seurat.annotationGraphs
		}

		annotation_graphs(input,file(assayParams.markersFile),assayname)
}
