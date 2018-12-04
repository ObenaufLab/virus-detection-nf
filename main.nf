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
        --inputDir        	Input directory of bam files.
        --output        	Output folder for centrifuge reports.

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
log.info " output directory         : ${params.output}"
log.info " ======================"
log.info ""
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

pairedEndRegex = params.inputDir + "/*_{1,2}.fq.gz"
SERegex = params.inputDir + "/*[!12].fq.gz"

pairFiles = Channel.fromFilePairs(pairedEndRegex)
singleFiles = Channel.fromFilePairs(SERegex, size: 1){ file -> file.baseName.replaceAll(/.fq/,"") }

singleFiles.mix(pairFiles)
.into { fastqPaVEChannel; fastqRefseqChannel; fastqENAChannel }
    
process centrifugePaVE {

	tag { lane }

    input:
    set val(lane), file(reads) from fastqPaVEChannel

    output:
    file ("*_PaVE_centrifuge_report.tsv") into PaVEChannel

    shell:
    
    def single = reads instanceof Path
    
	if (!single)
    
        '''
        
		centrifuge -x !{params.PaVEIndex} -q -p !{task.cpus} -1 !{reads[0]} -2 !{reads[1]} --report-file !{lane}_PaVE_centrifuge_report.tsv > /dev/null

        '''
    
    else 
    
        '''
        
		centrifuge -x !{params.PaVEIndex} -q -p !{task.cpus} -U !{reads} --report-file !{lane}_PaVE_centrifuge_report.tsv > /dev/null

        '''
}

process centrifugeRefSeq {

	tag { lane }

    input:
    set val(lane), file(reads) from fastqRefseqChannel

    output:
    file ("*_refseq_centrifuge_report.tsv") into refseqChannel

    shell:
    
    def single = reads instanceof Path
    
	if (!single)
    
        '''
        
		centrifuge -x !{params.refseqIndex} -q -p !{task.cpus} -1 !{reads[0]} -2 !{reads[1]} --report-file !{lane}_refseq_centrifuge_report.tsv > /dev/null

        '''
    
    else 
    
        '''
        
		centrifuge -x !{params.refseqIndex} -q -p !{task.cpus} -U !{reads} --report-file !{lane}_refseq_centrifuge_report.tsv > /dev/null

        '''
}

process centrifugeENA {

	tag { lane }

    input:
    set val(lane), file(reads) from fastqENAChannel

    output:
    file ("*_ENA_centrifuge_report.tsv") into ENAChannel

    shell:
    
    def single = reads instanceof Path
    
	if (!single)
    
        '''
        
		centrifuge -x !{params.ENAIndex} -q -p !{task.cpus} -1 !{reads[0]} -2 !{reads[1]} --report-file !{lane}_ENA_centrifuge_report.tsv > /dev/null

        '''
    
    else 
    
        '''
        
		centrifuge -x !{params.ENAIndex} -q -p !{task.cpus} -U !{reads} --report-file !{lane}_ENA_centrifuge_report.tsv > /dev/null

        '''
}