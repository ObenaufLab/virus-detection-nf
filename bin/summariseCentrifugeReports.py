#!/usr/bin/env python

from __future__ import print_function
import sys, os, glob

from argparse import ArgumentParser, RawDescriptionHelpFormatter

 # Info
usage = "Summarise centrifuge reports into a single report"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

parser.add_argument("-f", "--folder", type=str, required=True, dest="inputFolder", help="Folder with centrifuge report files")
parser.add_argument("-s", "--suffix", type=str, required=True, dest="suffix", help="Grouping suffix of centrifuge files - e.g. ENA_centrifuge_report.tsv. This will be trimmed from the basename of the file and used as sample name.")
parser.add_argument('-a', "--abundance", action='store_true', dest="abundance", help="Use abundance instead of unique reads in report.")
parser.add_argument('-t', "--taxID", action='store_true', dest="taxID", help="Use taxID instead of species name in report.")

args = parser.parse_args()

####################################################
# Get fastq gz filesizes
####################################################

taxIDToName = dict()
sampleToTaxa = dict()

reports = glob.glob(args.inputFolder + "/*" + args.suffix)

for report in reports:
#for x in xrange(0, 3):
    #report = reports[x]
    baseName = os.path.basename(report)[:-(len(args.suffix) + 1)]
    if not baseName in sampleToTaxa:
        sampleToTaxa[baseName] = dict()
    else :
        print("ERROR: Duplicate sample " + baseName + " in folder.", file = sys.stderr)
        sys.exit(-1)
    
    with open(report,'r') as f:
        next(f)
        for line in f:
            fields = line.rstrip().split("\t")
            name = fields[0]
            taxID = int(fields[1])
            reads = fields[5]
            abundance = fields[6]
            
            if not taxID in taxIDToName:
                taxIDToName[taxID] = name
                
            for sample in sampleToTaxa:
                if not taxID in sampleToTaxa[sample]:
                    sampleToTaxa[sample][taxID] = 0
     
            if args.abundance:
                sampleToTaxa[baseName][taxID] = abundance
            else :
                sampleToTaxa[baseName][taxID] = reads
           
sortedTaxa = taxIDToName.keys()
sortedTaxa.sort()

print("Sample",end="")

for taxon in sortedTaxa:
    if args.taxID:
        print("\t" + str(taxon),end="")
    else :
        print("\t" + taxIDToName[taxon],end="")
print()

for sample in sampleToTaxa:
    print(sample,end="")
    for taxon in sortedTaxa:
        if taxon in sampleToTaxa[sample]:
            print("\t" + str(sampleToTaxa[sample][taxon]),end="")
        else :
            print("\t0",end="")
    print()