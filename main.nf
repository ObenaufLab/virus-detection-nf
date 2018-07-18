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
log.info " Centrifuge index         : ${params.centrifugeIndex}"
log.info " Salmon index             : ${params.salmonIndex}"
log.info " bwa index                : ${params.bwaIndex}"
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
    
    input:
    set val(lane), file(bam) from rawBamFiles

    output:
    set val(lane), stdout, file("${lane}*.fq") into fastqFilesFromBamCentrifuge, fastqFilesFromBamSalmon, fastqFilesFromBwa

    shell:
    '''
    
    paired=`samtools view -c -f 1 !{bam}`
    
    if [ $paired -eq "0" ]; then
       	printf "False"
    	bamToFastq -i !{bam} -fq !{lane}.fq
    else
		printf "True"
		bamToFastq -i !{bam} -fq !{lane}_1.fq -fq2 !{lane}_2.fq
    fi
    
    '''
}

process centrifuge {

	tag { lane }
        
    input:
    set val(lane), val(paired), file(reads) from fastqFilesFromBamCentrifuge

    output:
    file ("*centrifuge_report.tsv") into centrifugeChannel

    shell:

    if( paired == 'True' )
        '''
	    centrifuge -x !{params.centrifugeIndex} -q -p !{task.cpus} -1 !{reads[0]} -2 !{reads[1]} --report-file !{lane}_centrifuge_report.tsv > /dev/null
	    '''
    else
        '''
	    centrifuge -x !{params.centrifugeIndex} -q -p !{task.cpus} -U !{reads} --report-file !{lane}_centrifuge_report.tsv > /dev/null
	    '''

}

process salmon {

	tag { lane }
        
    input:
    set val(lane), val(paired), file(reads) from fastqFilesFromBamSalmon

    output:
    file ("${lane}_salmon/quant.sf") into salmonChannel

    shell:

    if( paired == 'True' )
        '''
	    salmon quant -i !{params.salmonIndex} -l A -1 !{reads[0]} -2 !{reads[1]} -o !{lane}_salmon -p !{task.cpus}
	    '''
    else
        '''
	    salmon quant -i !{params.salmonIndex} -l A -r !{reads} -o !{lane}_salmon -p !{task.cpus}
	    '''

}

process bwa {

	tag { lane }
        
    input:
    set val(lane), val(paired), file(reads) from fastqFilesFromBwa

    output:
    set val(lane), file ("${lane}.bwa*") into bwaChannel

    shell:
    
    if (task.cpus > 8) {
    	bwaThreads = task.cpus / 4 * 3
    	sortThreads = task.cpus - bwaThreads
    } else {
    	bwaThreads = task.cpus
    	sortThreads = 1
    }
   
    '''
    bwa mem -M -R '@RG\\tID:!{lane}\\tPL:illumina\\tSM:!{lane}' -t !{bwaThreads} !{params.bwaIndex} !{reads} | \
    	samtools view -b - | samtools sort -@ !{sortThreads} -o !{lane}.bwa.bam
    	
    samtools index !{lane}.bwa.bam
    
    samtools idxstats !{lane}.bwa.bam > !{lane}.bwa.stats
    '''
}

process manta {

	tag { lane }
		     
    input:
    set val(lane), file(bwa) from bwaChannel
    
    output:
    file ("manta/results/variants/*") into outManta
    
    shell:
    '''
        
    shopt -s expand_aliases
    
    configManta.py --bam !{bwa[0]} \
    			   --referenceFasta !{params.bwaIndex} \
    			   --runDir manta    			   
    ${PWD}/manta/runWorkflow.py -m local -j !{task.cpus} -g !{task.memory.toGiga()}
	
    '''
}