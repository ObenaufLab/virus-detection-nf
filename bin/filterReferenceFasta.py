#!/usr/bin/env python

from __future__ import print_function
import sys, os

from Bio import SeqIO

from argparse import ArgumentParser, RawDescriptionHelpFormatter

 # Info
usage = "Filter duplicate fasta entries and uniquize them"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

parser.add_argument("-f", "--fasta", type=str, required=True, dest="fastaFile", help="Fasta file to design guides from")

args = parser.parse_args()

####################################################
# Finding duplicate entries
####################################################

print("Finding duplicate entries", file = sys.stderr)

chromosomes = dict()
duplicates = dict()

fasta_sequences = SeqIO.parse(open(args.fastaFile),'fasta')

for fasta in fasta_sequences:
    
    id, seq = fasta.id, str(fasta.seq)
    
    id = id.replace(",", "_")
        
    if id in chromosomes:
        duplicates[id] = 1
    else :
        chromosomes[id] = 0
        
print("Found " + str(len(duplicates)) + " duplicate entries.", file = sys.stderr)

####################################################
# Uniquizing entries
####################################################

print("Uniquizing entries.", file = sys.stderr)

fasta_sequences = SeqIO.parse(open(args.fastaFile),'fasta')

for fasta in fasta_sequences:
    
    id, seq = fasta.id, str(fasta.seq)
    
    id = id.replace(",", "_")
    
    if id in duplicates:
        print(">" + id + "_" + str(duplicates[id]))
        duplicates[id] += 1
    else :
        print(">" + id)
        
    print(seq)
    
print("Done.", file = sys.stderr)