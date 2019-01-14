#!/usr/bin/env python

from __future__ import print_function
import sys, os, re, gzip

from Bio import SeqIO

from argparse import ArgumentParser, RawDescriptionHelpFormatter

 # Info
usage = "Filter GenCode transcripts for cancer genes"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

parser.add_argument("-f", "--fasta", type=str, required=True, dest="fastaFile", help="Transcript fasta file")
parser.add_argument("-i", "--ids", type=str, required=True, dest="idFile", help="Ensembl transcript IDs to extract")

args = parser.parse_args()

####################################################
# Finding duplicate entries
####################################################

print("Reading transcript IDs to extract", file = sys.stderr)

txIds = list()

with open(args.idFile,'r') as f:
    for line in f:
        txIds.append(line.rstrip())

####################################################
# Extracting transcripts
####################################################

print("Extracting cancer gene transcripts.", file = sys.stderr)

fasta_sequences = SeqIO.parse(gzip.open(args.fastaFile),'fasta')

for fasta in fasta_sequences:
    
    id, seq = fasta.id, str(fasta.seq)
    
    txId = re.sub("\|.*","",id)
    
    if txId in txIds:
        print(">" + id + "\n" + seq)
    
print("Done.", file = sys.stderr)