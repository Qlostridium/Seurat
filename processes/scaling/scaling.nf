#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.folder.standAlone == true) ? "${params.folder.runDir}/processes": "${params.folder.runDir}/repos/Nf_Module_Seurat/processes"

process SEURAT__SCALING{
	//publishDir "${params.folder.outDir}/${samplename}", mode: 'symlink'
	input:
	tuple val(samplename), file(sobj)
	val assayType
	output:
	tuple val(samplename), file("${samplename}.SEURAT__SCALING_${assayType}.rds")
	script:
	def sampleParams = params.configParser(samplename, params.Seurat.scaling)
	def assayParams = params.configParser(assayType, sampleParams)
	"""
	Rscript ${scriptDir}/scaling/scaling.R --inputSeuratRds ${sobj} \
		--output "${samplename}.SEURAT__SCALING_${assayType}.rds" \
		--assay $assayType \
		${(!(assayParams.features == null)) ? '--features ' + assayParams.features: ''} \
		${(!(assayParams.varsToRegress == null)) ? '--varsToRegress ' + assayParams.varsToRegress: ''} \
		${(!(assayParams.splitBy == null)) ? '--splitBy ' + assayParams.splitBy: ''} \
		${(!(assayParams.modelUse == null)) ? '--modelUse ' + assayParams.modelUse: ''} \
		${(!(assayParams.useUmi == null)) ? '--useUmi ' + assayParams.useUmi: ''} \
		${(!(assayParams.doScale == null)) ? '--doScale ' + assayParams.doScale: ''} \
		${(!(assayParams.doCenter == null)) ? '--doCenter ' + assayParams.doCenter: ''} \
		${(!(assayParams.scaleMax == null)) ? '--scaleMax ' + assayParams.scaleMax: ''} \
		${(!(assayParams.blockSize == null)) ? '--blockSize ' + assayParams.blockSize: ''} \
		${(!(assayParams.minCellsToBlock == null)) ? '--minCellsToBlock ' + assayParams.minCellsToBlock: ''} \
		${(!(assayParams.verbose == null)) ? '--verbose ' + assayParams.verbose: ''}
	"""
}
