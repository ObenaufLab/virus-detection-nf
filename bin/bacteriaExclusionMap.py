#!/usr/bin/env python

from __future__ import print_function
import sys, os, re

from argparse import ArgumentParser, RawDescriptionHelpFormatter

 # Info
usage = "Create bacteria taxa exclusion map"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

parser.add_argument("-t", "--taxID", type=str, required=True, dest="taxIDMap", help="Map of sequence -> taxID of bacteria taxa")
parser.add_argument("-n", "--names", type=str, required=True, dest="nameMap", help="Map of taxID -> species names")
parser.add_argument("-d", "--nodes", type=str, required=True, dest="nodeMap", help="Map of taxID -> ")

args = parser.parse_args()

####################################################
# Get tax IDs
####################################################

taxa = dict()

with open(args.taxIDMap, 'r') as f:
    for line in f:
        sepLine = re.sub(r'\s+','\t',line.rstrip())
        seq, taxon = sepLine.split("\t")
        taxa[taxon] = ""
        
####################################################
# Fetch all names
####################################################

names = dict()
        
with open(args.nameMap, 'r') as f:
    for line in f:
        fields = line.rstrip().split("|")
        
        id = fields[0]
        id = re.sub("^\t+","",id)
        id = re.sub("\t+$","",id)
        
        name = fields[1]
        name = re.sub("^\t+","",name)
        name = re.sub("\t+$","",name)
        
        classification = fields[3]
        classification = re.sub("^\t+","",classification)
        classification = re.sub("\t+$","",classification)
        
        if not id in names:
            names[id] = name
        if classification == "scientific name":
            names[id] = name
            
####################################################
# Recursively get all higher order taxa
####################################################

fullTaxonomy = dict()
curLength = -1
iter = 0

while curLength != len(fullTaxonomy):
    iter += 1
    print("Taxonomy iteration: " + str(iter), file=sys.stderr)
    curLength = len(fullTaxonomy)
    with open(args.nodeMap, 'r') as f:
        for line in f:
            fields = line.rstrip().split("|")
            
            leaf = fields[0]
            leaf = re.sub("^\t+","",leaf)
            leaf = re.sub("\t+$","",leaf)
            
            ancestor = fields[1]
            ancestor = re.sub("^\t+","",ancestor)
            ancestor = re.sub("\t+$","",ancestor)
            
            if (leaf in taxa) :
                
                #taxa[leaf] = ""
                taxa[ancestor] = ""
            
                if not leaf in fullTaxonomy:
                    fullTaxonomy[leaf] = names[leaf]
                if not ancestor in fullTaxonomy:
                    fullTaxonomy[ancestor] = names[ancestor]

for taxon in fullTaxonomy:
    print(taxon + "\t" + fullTaxonomy[taxon])
