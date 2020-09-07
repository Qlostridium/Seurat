#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (!(params.folder.standAlone == null) && params.folder.standAlone == true) ? "${params.global.rundir}/processes": "${params.global.rundir}/src/Seurat/processes"

process SEURAT__HTO_DEMULTIPLEXNG {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	cache 'lenient'
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink', pattern : "${samplename}_logQC.txt"
	container params.Seurat.container
  input:
	tuple val(samplename), file(seuratobj)
	file(freemuxfile)
	file(qcfile)
  output:
	tuple val(samplename), file("${samplename}.SEURAT__HTO_DEMULTIPLEXING.rds")
	file("${samplename}_logQC.txt")
	file("${samplename}_confusionmatrix.xlsx")
  script:
	def realqcfile = qcfile.name != 'NO_FILE' ? "--qcFile $qcfile" : ''
	"""
	Rscript ${scriptDir}/HTODemultiplexing/HTODemultiplexing.R --seuratObj ${seuratobj} \
	--output "${samplename}.SEURAT__HTO_DEMULTIPLEXING.rds" \
	--scriptFunctions ${scriptDir}/utils/script_functions_COVID.R \
	--freemuxFile $freemuxfile \
	$realqcfile
	"""
}
