#!/usr/bin/env python

from __future__ import print_function
import sys, os, glob

from argparse import ArgumentParser, RawDescriptionHelpFormatter

 # Info
usage = "Compare sorted aln bam vs gzipped fastq compression"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

parser.add_argument("-f", "--fastq", type=str, required=True, dest="fastqFolder", help="Folder with fastq files")
parser.add_argument("-b", "--bam", type=str, required=True, dest="bamFolder", help="Folder with corresponding bam files")

args = parser.parse_args()

####################################################
# Get fastq gz filesizes
####################################################

fastqs = glob.glob(args.fastqFolder + "/*fq.gz")

samples = dict()

for file in fastqs:
    baseName = os.path.basename(file)[:-8]
    fileSize = os.path.getsize(file) >> 20
    
    if not baseName in samples:
        samples[baseName] = 0
    samples[baseName] += fileSize
    
####################################################
# Compare to BAMs
####################################################

print("File\tBAM_size\tFastqGz_size")

bams = glob.glob(args.bamFolder + "/*/*bam")

for file in bams:
    baseName = os.path.basename(file)[:-4]
    fileSize = os.path.getsize(file) >> 20
    
    if baseName in samples:
        print(baseName + "\t" + str(fileSize) + "\t" + str(samples[baseName]))