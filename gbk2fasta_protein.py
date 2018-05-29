#!/usr/bin/env python


# Created:               28-05-2018                  gabuali
# Last Modified:         Mon May 28 16:43:06 2018    gabuali


# Extract amino acid sequence from GBK file and save in FASTA format.
# Header contains gene locus_tag, protein ID, and functional characterization.

from sys import argv
from Bio import SeqIO

gbk_in = argv[1]
faa_out = gbk_in.replace( "genomic.gbff", "protein.faa" )
fh_in  = open(gbk_in, "r")
fh_out = open(faa_out, "w")


# get all sequence records for the genbank file
recs = [ rec for rec in SeqIO.parse( argv[1], "genbank" ) ]

# iterate over each sequence record (contig in draft genome) and extract info,
# skipping features source, gene, operon
for rec in recs:

# NOTE that Main features, eg  accession No, etc. are NOT provided in Prokka annotations
    print( len(recs) )                # number of extracted contigs or sequence records
    organism = rec.description          # DEFINITION field

    if rec.features:
        feats = rec.features

        for feat in feats:
            if feat.type == "CDS":
                if not 'pseudo' in feat.qualifiers.keys():
                    # feature type
                    f_type = feat.type # feature type (e.g. CDS, rRNA, etc.)

                    # not all features have all qualifiers so
                    # if they're absent add 'NA'
                    if 'locus_tag' in feat.qualifiers.keys():
                        locus = feat.qualifiers['locus_tag'][0] # locus_tag
                    else:
                        locus = 'NA'

                    if 'note' in feat.qualifiers.keys():
                        note = feat.qualifiers['note'][0]
                    else:
                        note = 'NA'

                    if 'product' in feat.qualifiers.keys():
                        product = feat.qualifiers['product'][0] # functional characterization
                    else:
                        product = 'NA'

                    if 'protein_id' in feat.qualifiers.keys():
                        protein_id = feat.qualifiers['protein_id'][0]
                    else:
                        protein_id = 'NA'

                    if 'translation' in feat.qualifiers.keys():
                        translation = feat.qualifiers['translation'][0] # aa sequence
                    else:
                        translation = 'NA'

                    fh_out.write( f">{locus}__{protein_id}__{product}\n{translation}\n" )


fh_in.close()
fh_out.close()
