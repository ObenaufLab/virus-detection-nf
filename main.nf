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
    nextflow run obenauflab/virus-detection-nf -r manta

    Options:
        --inputDir        	Input directory of fastq files.
        --outputDir        	Output folder for SV vcf files.
        --index             Index to use, one of PaVE|RefSeq|ENA

    Profiles:
        standard            local execution
        ii2                 SLURM execution with singularity on IMPIMBA2
        aws                 SLURM execution with singularity on IMPIMBA2

    Docker:
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
log.info " output directory         : ${params.outputDir}"
log.info " ======================"
log.info ""

if (params.index == "PaVE") {

  centrifugeIndex = Channel
	   .fromPath(params.centrifugePaVEIndex)
	   .ifEmpty { exit 1, "Centrifuge index not found: ${params.centrifugePaVEIndex}" }

   bwaIndex = Channel
 	   .fromPath(params.BWAPaVEIndex)
 	   .ifEmpty { exit 1, "BWA index not found: ${params.BWAPaVEIndex}" }

   mantaIndex = Channel
       .fromPath(params.BWAPaVEIndex)
       .ifEmpty { exit 1, "Manta index not found: ${params.BWAPaVEIndex}" }

} else if (params.index == "RefSeq") {

  centrifugeIndex = Channel
    .fromPath(params.centrifugeRefSeqIndex)
    .ifEmpty { exit 1, "Centrifuge index not found: ${params.centrifugeRefSeqIndex}" }

   bwaIndex = Channel
      .fromPath(params.BWARefSeqIndex)
      .ifEmpty { exit 1, "BWA index not found: ${params.BWARefSeqIndex}" }

   mantaIndex = Channel
       .fromPath(params.BWARefSeqIndex)
       .ifEmpty { exit 1, "Manta index not found: ${params.BWARefSeqIndex}" }

} else if (params.index == "ENA") {

  centrifugeIndex = Channel
    .fromPath(params.centrifugeENAIndex)
    .ifEmpty { exit 1, "Centrifuge index not found: ${params.centrifugeENAIndex}" }

   bwaIndex = Channel
      .fromPath(params.BWAENAIndex)
      .ifEmpty { exit 1, "BWA index not found: ${params.BWAENAIndex}" }

   mantaIndex = Channel
       .fromPath(params.BWAENAIndex)
       .ifEmpty { exit 1, "Manta index not found: ${params.BWAENAIndex}" }

} else {

   log.info"""
   Choose one of PaVE | RefSeq | ENA as Index!

   """.stripIndent()
   helpMessage()
   exit 0

}

pairedEndRegex = params.inputDir + "/*_{1,2}.fq.gz"
SERegex = params.inputDir + "/*[!12].fq.gz"

pairFiles = Channel.fromFilePairs(pairedEndRegex)
singleFiles = Channel.fromFilePairs(SERegex, size: 1){ file -> file.baseName.replaceAll(/.fq/,"") }

singleFiles.mix(pairFiles)
.set { fastqChannel }

mantaConfigFile = file(params.mantaConfig)

/*process centrifugeMatchExtraction {

	tag { lane }

    input:
    set val(lane), file(reads) from fastqChannel
    file index from centrifugeIndex.first()


    output:
    set val(lane), file ("reads*fq") into readSubsetChannel

    shell:

    def single = reads instanceof Path

	if (!single)

        '''
        centrifuge -x !{index}/centrifuge_index -q -p !{task.cpus} -1 !{reads[0]} -2 !{reads[1]} | grep 20000109 | cut -f 1 | sort | uniq > readIDs
        seqtk subseq !{reads[0]} readIDs > reads1.fq
        seqtk subseq !{reads[1]} readIDs > reads2.fq

        '''

    else

        '''

        centrifuge -x !{index}/centrifuge_index -q -p !{task.cpus} -U !{reads} | grep 20000109 | cut -f 1 | sort | uniq > readIDs
        seqtk subseq !{reads} readIDs > reads.fq

        '''
}*/

process bwa {

	tag { lane }

    input:
    // set val(lane), file(reads) from readSubsetChannel
    set val(lane), file(reads) from fastqChannel
    file index from bwaIndex.first()

    output:
    set val(lane), file ("${lane}.bwa*") into bwaChannel

    shell:

    if (task.cpus > 4) {
    	bwaThreads = task.cpus - 2
    	sortThreads = 2
    } else {
    	bwaThreads = task.cpus
    	sortThreads = 1
    }

    '''
    bwa mem -M -R '@RG\\tID:!{lane}\\tPL:illumina\\tSM:!{lane}' -t !{bwaThreads} !{index}/bwa_index.fa !{reads} | \
    	samtools view -b - | samtools sort -@ !{sortThreads} -o !{lane}.bwa.bam

    samtools index !{lane}.bwa.bam
    '''
}

process manta {

	tag { lane }

    input:
    set val(lane), file(bwa) from bwaChannel
    file mantaConfigFile
    file index from mantaIndex.first()

    output:
    file ("manta/results/variants/*") into outManta

    shell:

    def configArg = mantaConfigFile.name != 'NO CONFIG' ? "--config ${params.mantaConfig}" : ''
    '''

    shopt -s expand_aliases

    configManta.py --bam !{bwa[0]} \
    			   --referenceFasta !{index}/bwa_index.fa \
    			   --runDir manta \
             --rna \
             --generateEvidenceBam !{configArg}

    ${PWD}/manta/runWorkflow.py -m local -j !{task.cpus} -g !{task.memory.toGiga()}

    mkdir -p !{lane}
    mv manta/results/variants/candidateSmallIndels.vcf.gz manta/results/variants/!{lane}_candidateSmallIndels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi manta/results/variants/!{lane}_candidateSmallIndels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz manta/results/variants/!{lane}_candidateSV.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi manta/results/variants/!{lane}_candidateSV.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz manta/results/variants/!{lane}_diploidSV.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi manta/results/variants/!{lane}_diploidSV.vcf.gz.tbi

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
