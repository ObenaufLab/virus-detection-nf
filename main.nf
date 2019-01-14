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
    nextflow run obenauflab/virus-detection-nf -r salmon

    Options:
        --inputDir        	Input directory of fastq files.
        --outputDir        	Output folder for salmon quantification files.

    Profiles:
        standard            local execution
        ii2                 SLURM execution with singularity on IMPIMBA2
        aws                 SLURM execution with singularity on IMPIMBA2

    Docker:
    obenauflab/salmon:latest

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
log.info " output directory         : ${params.output}"
log.info " ======================"
log.info ""

log.info ""
log.info " parameters "
log.info " ======================"
log.info " input directory          : ${params.inputDir}"
log.info " output directory         : ${params.outputDir}"
log.info " ======================"
log.info ""

pairedEndRegex = params.inputDir + "/*_{1,2}.fq.gz"
SERegex = params.inputDir + "/*[!12].fq.gz"

pairFiles = Channel.fromFilePairs(pairedEndRegex)
singleFiles = Channel.fromFilePairs(SERegex, size: 1){ file -> file.baseName.replaceAll(/.fq/,"") }

singleFiles.mix(pairFiles)
.set { fastqChannel }

indexChannel = Channel
	.fromPath(params.salmonIndex)
	.ifEmpty { exit 1, "Salmon index not found: ${params.salmonIndex}" }

process salmon {

	tag { lane }

    input:
    set val(lane), file(reads) from fastqChannel
    file index from indexChannel.first()

    output:
    file ("${lane}_salmon/quant.sf") into salmonChannel
    file ("${lane}_pseudo.bam") into pseudoBamChannel

    shell:

    def single = reads instanceof Path

    if (!single)

      '''
      salmon quant -i !{index} -l A -1 !{reads[0]} -2 !{reads[1]} -z -o !{lane}_salmon -p !{task.cpus} | samtools view -Sb -F 256 - > !{lane}_pseudo.bam
	    '''
    else
      '''
      salmon quant -i !{index} -l A -r !{reads} -z -o !{lane}_salmon -p !{task.cpus} | samtools view -Sb -F 256 - > !{lane}_pseudo.bam
	    '''

}

workflow.onComplete {
	RED='\033[0;31m'
    GREEN='\033[0;32m'
    NC='\033[0m'

    log.info "\nobenauflab/virus-detection has finished."
    log.info "Status:   " + (workflow.success ? "${GREEN}SUCCESS${NC}" : "${RED}ERROR${NC}")
    log.info "Time:     ${workflow.complete}"
    log.info "Duration: ${workflow.duration}\n"
}
