#!/usr/bin/python
import pysam

import argparse
parser = argparse.ArgumentParser(description='Small script to add custom IDs to a VCF file')
parser.add_argument("--i", default=1, help="This is the input vcf file")
parser.add_argument("--o", default=1, help="This is the output vcf file")

args = parser.parse_args()
infile = args.i
outfile = args.o

#print(infile)
#print(outfile)

# read the input file
myvcf = pysam.VariantFile(infile, "r")

# Add the field to header. Say its a string and can take any values.
myvcf.header.formats.add("MERGEID", ".", "String", "ID to merge MAF and VCF after conversion")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile(outfile, 'w', header=myvcf.header)

# iterate over all variants and add the custom defined header
with open(outfile, "a") as out:
    out.write(str(myvcf.header))
    for variant in myvcf:
        for sample in variant.samples:
            newEntry=variant.contig + "-" + str(variant.pos) + "-" + variant.ref
            variant.samples[sample]['MERGEID']=newEntry
            out.write(str(variant))

#variant.contig + ":" + str(variant.pos) + ":" + variant.ref + ":" + variant.alts
