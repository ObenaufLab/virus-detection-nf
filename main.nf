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
    nextflow run obenauflab/virus-detection-nf -r trinity

    Options:
        --inputDir        	Input directory of bam files.
        --outputDir        	Output folder for Trinity.

    Profiles:
        standard            local execution
        ii2                 SLURM execution with singularity on IMPIMBA2
        aws                 AWS batch execution

    Docker:
    quay.io/biocontainers/samtools:1.9--h8ee4bcc_1
    trinityrnaseq/trinityrnaseq:2.8.4

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

    	 samtools fastq -@ !{task.cpus} -s !{lane}.fq !{bam}

    else

       printf "True"

		   samtools collate -f -O -u -@ !{task.cpus} !{bam} | samtools fastq -1 !{lane}_1.fq -2 !{lane}_2.fq -N -@ !{task.cpus} -

    fi

    '''
}

process trinity {

    label 'trinity'

    tag { lane }

    publishDir path: "${params.outputDir}/trinity/", mode: 'copy',
               overwrite: 'true', pattern: "*trinity.fa"

    input:
    set val(lane), val(paired), file(reads) from rawReads

    output:
    set val(lane), file("*trinity.fa") into trinityTransdecoder, trinityBlastx, trinityTrinotate
    set val(lane), file("*trinity.fa.gene_trans_map") into outTrinityTransMap

    shell:
    if( paired == 'True' )
      '''
      Trinity --seqType fq --SS_lib_type RF \
               --max_memory !{task.memory.toGiga()}G \
               --left !{reads[0]} \
               --right !{reads[1]} \
               --CPU !{task.cpus} --output trinity

      mv trinity/Trinity.fasta !{lane}_trinity.fa

      /usr/local/bin/trinityrnaseq/util/support_scripts/get_Trinity_gene_to_trans_map.pl !{lane}_trinity.fa > !{lane}_trinity.fa.gene_trans_map
	    '''
    else
      '''
      Trinity --seqType fq --SS_lib_type RF \
               --max_memory !{task.memory.toGiga()}G \
               --single !{reads} \
               --CPU !{task.cpus} --output trinity

      mv trinity/Trinity.fasta !{lane}_trinity.fa

      /usr/local/bin/trinityrnaseq/util/support_scripts/get_Trinity_gene_to_trans_map.pl !{lane}_trinity.fa > !{lane}_trinity.fa.gene_trans_map
	    '''
}

process transdecoder {

    tag { lane }

    input:
    set val(lane), file(transcripts) from trinityTransdecoder

    output:
    set val(lane), file("*trinity.fa.transdecoder.pep") into transdecoderBlastp, transdecoderHmmscan, transdecoderTrinotate

    shell:
    '''
    TransDecoder.LongOrfs -t !{transcripts}
    TransDecoder.Predict -t !{transcripts}
    '''
}

process blastx {

    label 'trinotate'

    tag { lane }

    input:
    set val(lane), file(transcripts) from trinityBlastx

    output:
    set val(lane), file("*blastx.outfmt6") into outBlastx

    shell:
    '''
    blastx -query !{transcripts} -db !{params.blastdb} -num_threads !{task.cpus} -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > !{lane}_blastx.outfmt6
	  '''
}

process blastp {

    label 'trinotate'

    tag { lane }

    input:
    set val(lane), file(proteins) from transdecoderBlastp

    output:
    set val(lane), file("*blastp.outfmt6") into outBlastp

    shell:
    '''
    blastp -query !{proteins} -db !{params.blastdb} -num_threads !{task.cpus} -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > !{lane}_blastp.outfmt6
	  '''
}

process hmmscan {

    label 'trinotate'

    tag { lane }

    input:
    set val(lane), file(proteins) from transdecoderHmmscan

    output:
    set val(lane), file("*TrinotatePFAM.out") into outHmmScan

    shell:
    '''
    hmmscan --cpu !{task.cpus} --domtblout !{lane}_TrinotatePFAM.out !{params.pfamdb} !{proteins} > pfam.log
	  '''
}

process trinotatedb {

    label 'trinotate'

    tag { lane }

    publishDir path: "${params.outputDir}/trinotate/", mode: 'copy',
               overwrite: 'true', pattern: "*trinotate_annotation_report.txt"

    input:
    set val(lane), file(transcripts) from trinityTrinotate
    set val(lane), file(transMap) from outTrinityTransMap
    set val(lane), file(proteins) from transdecoderTrinotate
    set val(lane), file(blastx) from outBlastx
    set val(lane), file(blastp) from outBlastp
    set val(lane), file(hmmscan) from outHmmScan
    file(sqlite) from sqlite.collect()

    output:
    set val(lane), file("*trinotate_annotation_report.txt") into outTrinotateDb

    shell:
    '''
    cp !{params.sqlite} Trinotate.sqlite
    Trinotate Trinotate.sqlite init --gene_trans_map !{transMap} --transcript_fasta !{transcripts} --transdecoder_pep !{proteins}
    Trinotate Trinotate.sqlite LOAD_swissprot_blastp !{blastp}
    Trinotate Trinotate.sqlite LOAD_swissprot_blastx !{blastx}
    Trinotate Trinotate.sqlite LOAD_pfam !{hmmscan}
    Trinotate Trinotate.sqlite report > !{lane}_trinotate_annotation_report.txt
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
