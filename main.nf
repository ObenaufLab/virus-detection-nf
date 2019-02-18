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
    nextflow run obenauflab/virus-detection-nf -r STARfusion

    Options:
        --inputDir        	Input directory of bam files.
        --outputDir        	Output folder for STAR-fusion results.

    Profiles:
        standard            local execution
        ii2                 SLURM execution with singularity on IMPIMBA2
        aws                 AWS batch execution

    Docker:
    quay.io/biocontainers/samtools:1.9--h8ee4bcc_1
    trinityctat/ctatfusion:1.5.0

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

Channel
    .fromPath( "${params.inputDir}/*.bam" )
    .map { file -> tuple( file.baseName, file ) }
    .set { rawBamFiles }

process bamToFastq {

    tag { lane }

    input:
    set val(lane), file(bam) from rawBamFiles

    output:
    set val(lane), stdout, file("${lane}*.fq") into rawReads

    shell:
    '''

    paired=`samtools view -c -f 1 !{bam}`

    if [ $paired -eq "0" ]; then

       printf "False"

    	 samtools fastq -@ !{task.cpus} -0 !{lane}.fq !{bam}

    else

       printf "True"

		   samtools collate -f -O -u -@ !{task.cpus} !{bam} | samtools fastq -1 !{lane}_1.fq -2 !{lane}_2.fq -N -@ !{task.cpus} -

    fi

    '''
}

process starfusion {

    errorStrategy 'ignore'

    tag { lane }

    input:
    set val(lane), val(paired), file(reads) from rawReads

    output:
    set val(lane), file("*_STARFusion"), file("*finspector.txt"), file("*trinity_fusions.fa"), file("*trinity_fusions.bed.gz") into outTrinity

    shell:
    if( paired == 'True' )
      '''
      /usr/local/src/STAR-Fusion/STAR-Fusion \
               --genome_lib_dir !{params.STARFusionIndex} \
               --left_fq !{reads[0]} \
               --right_fq !{reads[1]} \
               --output_dir star_fusion_outdir \
               --FusionInspector validate \
               --denovo_reconstruct \
               --examine_coding_effect \
               --CPU !{task.cpus}

      mkdir -p !{lane}_STARFusion
      cp star_fusion_outdir/FusionInspector-validate/finspector.fusion_predictions.final.abridged.FFPM.annotated.coding_effect !{lane}_STARFusion/!{lane}_finspector.txt

       mv star_fusion_outdir/FusionInspector-validate/finspector.fusion_predictions.final.abridged.FFPM.annotated.coding_effect !{lane}_finspector.txt
       mv star_fusion_outdir/FusionInspector-validate/finspector.gmap_trinity_GG.fusions.fasta !{lane}_trinity_fusions.fa
       mv star_fusion_outdir/FusionInspector-validate/finspector.gmap_trinity_GG.fusions.gff3.bed.sorted.bed.gz !{lane}_trinity_fusions.bed.gz
       mv star_fusion_outdir/FusionInspector-validate/finspector.fa !{lane}_finspector.fa
       mv star_fusion_outdir/FusionInspector-validate/finspector.gtf !{lane}_finspector.gtf
       mv star_fusion_outdir/FusionInspector-validate/finspector.junction_reads.bam !{lane}_finspector.junction_reads.bam
       mv star_fusion_outdir/FusionInspector-validate/finspector.junction_reads.bam.bai !{lane}_finspector.junction_reads.bam.bai
       mv star_fusion_outdir/FusionInspector-validate/finspector.spanning_reads.bam !{lane}_finspector.spanning_reads.bam
       mv star_fusion_outdir/FusionInspector-validate/finspector.spanning_reads.bam.bai !{lane}_finspector.spanning_reads.bam.bai
	    '''
    else
      '''
      /usr/local/src/STAR-Fusion/STAR-Fusion \
               --genome_lib_dir !{params.STARFusionIndex} \
               --left_fq !{reads} \
               --output_dir star_fusion_outdir \
               --FusionInspector validate \
               --denovo_reconstruct \
               --examine_coding_effect \
               --CPU !{task.cpus}



       mv star_fusion_outdir/FusionInspector-validate/finspector.fusion_predictions.final.abridged.FFPM.annotated.coding_effect !{lane}_finspector.txt
       mv star_fusion_outdir/FusionInspector-validate/finspector.gmap_trinity_GG.fusions.fasta !{lane}_trinity_fusions.fa
       mv star_fusion_outdir/FusionInspector-validate/finspector.gmap_trinity_GG.fusions.gff3.bed.sorted.bed.gz !{lane}_trinity_fusions.bed.gz
       mv star_fusion_outdir/FusionInspector-validate/finspector.fa !{lane}_finspector.fa
       mv star_fusion_outdir/FusionInspector-validate/finspector.gtf !{lane}_finspector.gtf
       mv star_fusion_outdir/FusionInspector-validate/finspector.junction_reads.bam !{lane}_finspector.junction_reads.bam
       mv star_fusion_outdir/FusionInspector-validate/finspector.junction_reads.bam.bai !{lane}_finspector.junction_reads.bam.bai
       mv star_fusion_outdir/FusionInspector-validate/finspector.spanning_reads.bam !{lane}_finspector.spanning_reads.bam
       mv star_fusion_outdir/FusionInspector-validate/finspector.spanning_reads.bam.bai !{lane}_finspector.spanning_reads.bam.bai

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
