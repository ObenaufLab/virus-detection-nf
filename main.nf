#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2018 Tobias Neumann
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

def helpMessage() {
    log.info"""
    ================================================================
    virus-detection-nf
    ================================================================
    DESCRIPTION
    Usage:
    nextflow run obenauflab/virus-detection-nf
    Options:
        --inputDir        	Input directory of bam files.

    Profiles:
        standard            local execution
        sge			        SGE execution with singularity on IMPIMBA1
        ii2                 SLURM execution with singularity on IMPIMBA2
        
    Docker:
    combinelab/salmon:0.8.1
    obenauflab/virusintegration:latest
    
    Author:
    Tobias Neumann (tobias.neumann@imp.ac.at)
    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

log.info ""
log.info " parameters "
log.info " ======================"
log.info " input directory          : ${params.inputDir}"
log.info " ======================"
log.info ""
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

Channel
    .fromPath( "${params.inputDir}/*.bam" )
    .map { file -> tuple( file.baseName, file ) }
    .set { rawBamFiles }
    
process bamToFastq {

    tag { lane }
    
    container = 'docker://obenauflab/virusintegration:latest'

    input:
    set val(lane), file(bam) from rawBamFiles

    output:
    set val(lane), file("${lane}.fastq.gz") into fastqFilesFromBam

    shell:
    '''
    
    PAIRED=`samtools view -c -f 1 !{bam}`
    
    if $PAIRED == 0; then
    	echo "MUAHAHAHA"
    else
    	echo "NOOOO"
    fi
    
    #samtools fastq ${bam} | gzip -c > ${lane}.fastq.gz
    '''
}