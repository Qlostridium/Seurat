#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process annotation_graphs {
	publishDir "${params.folder.outDir}/${params.sampleName}", mode: 'symlink', pattern: "Plots/**"
	publishDir "${params.folder.outDir}/${params.sampleName}", mode: 'symlink', pattern: "allClusters_${samplename}.xlsx"
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

		//if(params.Seurat.annotationGraphs.containsKey(input[0]))
		//	assayParams = sampleParams."${input[0]}"
		//else
		//	assayParams = sampleParams

		if(params.Seurat.annotationGraphs.containsKey(assayname)){
			assayParams = params.Seurat.annotationGraphs."${assayname}"
		} else {
			assayParams = params.Seurat.annotationGraphs
		}

		annotation_graphs(input,file(assayParams.markersFile),assayname)
}
