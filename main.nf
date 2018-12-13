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
    virus-detection-nf centrifuge
    ================================================================
    DESCRIPTION
    Usage:
    nextflow run obenauflab/virus-detection-nf -r centrifuge
    Options:
        --inputDir        	Input directory of fastq files.
        --outputDir        	Output folder for centrifuge reports.

    Profiles:
        standard            local execution
        slurm			    SLURM execution with singularity on IMPIMBA2
        awsbatch            AWS batch execution
        
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

pairedEndRegex = params.inputDir + "/*_{1,2}.fq.gz"
SERegex = params.inputDir + "/*[!12].fq.gz"

pairFiles = Channel.fromFilePairs(pairedEndRegex)
singleFiles = Channel.fromFilePairs(SERegex, size: 1){ file -> file.baseName.replaceAll(/.fq/,"") }

singleFiles.mix(pairFiles)
.into { fastqPaVEChannel; fastqRefseqChannel; fastqENAChannel }

PaVEIndex = Channel
	.fromPath(params.PaVEIndex)
	.ifEmpty { exit 1, "PaVE index not found: ${params.PaVEIndex}" }
	
RefseqIndex = Channel
	.fromPath(params.refseqIndex)
	.ifEmpty { exit 1, "Refseq index not found: ${params.refseqIndex}" }
	
ENAIndex = Channel
	.fromPath(params.ENAIndex)
	.ifEmpty { exit 1, "ENA index not found: ${params.ENAIndex}" }
    
process centrifugePaVE {

	tag { lane }

    input:
    set val(lane), file(reads) from fastqPaVEChannel
    file index from PaVEIndex.first()
    

    output:
    file ("*_PaVE_centrifuge_report.tsv") into PaVEChannel
    file ("*_PaVE_kreport.out") into PaVEKReportChannel

    shell:
    
    def single = reads instanceof Path
    
	if (!single)
    
        '''
        
		centrifuge -x !{index}/obenauf_PaVE_index -q -p !{task.cpus} -1 !{reads[0]} -2 !{reads[1]} --report-file !{lane}_PaVE_centrifuge_report.tsv > readwise.txt
		centrifuge-kreport -x !{index}/obenauf_PaVE_index readwise.txt > !{lane}_PaVE_kreport.out

        '''
    
    else 
    
        '''
        
		centrifuge -x !{index}/obenauf_PaVE_index -q -p !{task.cpus} -U !{reads} --report-file !{lane}_PaVE_centrifuge_report.tsv > readwise.txt
		centrifuge-kreport -x !{index}/obenauf_PaVE_index readwise.txt > !{lane}_PaVE_kreport.out

        '''
}

process centrifugeRefSeq {

	tag { lane }

    input:
    set val(lane), file(reads) from fastqRefseqChannel
    file index from RefseqIndex.first()

    output:
    file ("*_refseq_centrifuge_report.tsv") into refseqChannel
    file ("*_refseq_kreport.out") into refseqKReportChannel

    shell:
    
    def single = reads instanceof Path
    
	if (!single)
    
        '''
        
		centrifuge -x !{index}/obenauf_plain_index -q -p !{task.cpus} -1 !{reads[0]} -2 !{reads[1]} --report-file !{lane}_refseq_centrifuge_report.tsv > readwise.txt
		centrifuge-kreport -x !{index}/obenauf_plain_index readwise.txt > !{lane}_refseq_kreport.out

        '''
    
    else 
    
        '''
        
		centrifuge -x !{index}/obenauf_plain_index -q -p !{task.cpus} -U !{reads} --report-file !{lane}_refseq_centrifuge_report.tsv > readwise.txt
		centrifuge-kreport -x !{index}/obenauf_plain_index readwise.txt > !{lane}_refseq_kreport.out

        '''
}

process centrifugeENA {

	tag { lane }

    input:
    set val(lane), file(reads) from fastqENAChannel
    file index from ENAIndex.first()

    output:
    file ("*_ENA_centrifuge_report.tsv") into ENAChannel
    file ("*_ENA_kreport.out") into ENAKReportChannel

    shell:
    
    def single = reads instanceof Path
    
	if (!single)
    
        '''
        
		centrifuge -x !{index}/obenauf_ENA_index -q -p !{task.cpus} -1 !{reads[0]} -2 !{reads[1]} --report-file !{lane}_ENA_centrifuge_report.tsv > readwise.txt
		centrifuge-kreport -x !{index}/obenauf_ENA_index readwise.txt > !{lane}_ENA_kreport.out

        '''
    
    else 
    
        '''
        
		centrifuge -x !{index}/obenauf_ENA_index -q -p !{task.cpus} -U !{reads} --report-file !{lane}_ENA_centrifuge_report.tsv > readwise.txt
		centrifuge-kreport -x !{index}/obenauf_ENA_index readwise.txt > !{lane}_ENA_kreport.out

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