#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include {run_HTO } from'./workflows/HTO.nf' params(params)
include {run_RNA } from'./workflows/RNA.nf' params(params)
include {run_ADT } from'./workflows/ADT.nf' params(params)

workflow HTO {
	input = Channel.fromPath(params.Seurat.seuratObjBuilder.inputFile)
					.map{ file -> tuple(params.sampleName,file)}
	run_HTO(input)
}

workflow RNA {
	main:
	if(params.Seurat.seuratObjBuilder.inputFile == null){
		seuratInput = Channel.fromPath(params.inputFile)
						.map{ file -> tuple(params.sampleName,file)}
	} else {
		input = Channel.fromPath(params.Seurat.seuratObjBuilder.inputFile)
						.map{ file -> tuple(params.sampleName,file)}
		seuratInput = SEURAT__SEURAT_OBJECT_BUILDER(input)
	}
	run_RNA(seuratInput)
	emit:
	run_RNA.out
}

workflow ADT {
	input = Channel.fromPath(params.inputFile)
					.map{ file -> tuple(params.sampleName,file)}
	run_ADT(input)
}

workflow hashtags_rnaseq {
	input = Channel.fromPath(params.Seurat.seuratObjBuilder.inputFile)
					.map{ file -> tuple(params.sampleName,file)}
	run_HTO(input)
	run_RNA(run_HTO.out)
}

workflow hashtags_citeseq {
	input = Channel.fromPath(params.Seurat.seuratObjBuilder.inputFile)
					.map{ file -> tuple(params.sampleName,file)}
	run_HTO(input)
	run_RNA(run_HTO.out)
	run_ADT(run_RNA.out)
}

workflow citeseq {
	RNA()
	run_ADT(RNA.out)
}
