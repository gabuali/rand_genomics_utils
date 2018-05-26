#!/usr/bin/env python


# Created:               17-05-2018                  gabuali
# Last Modified:         Fri May 18 07:36:27 2018    gabuali

# Parse genbank record and convert to tab-delimited file.
# Tested on bacterial genbank records, both complete and draft.
# Feature qualifiers specific to eukaryotic records are not captured here but can 
# easily be added.

import argparse
import os
from Bio import SeqIO


# The if condition will execute argparse only if script (module) is run directly,
# and not when imported into another program.
if __name__ == '__main__':
    # define input args and parse them from the command line
    parser = argparse.ArgumentParser( description = "Convert genbank record into tab-separated table." )

    # input and output files as required named arguments
    requiredNamed = parser.add_argument_group( 'required named arguments' )
    requiredNamed.add_argument( '--infile', '-i', metavar = 'GENBANK', type = argparse.FileType( 'r', encoding = 'UTF-8' ), required = True, help = "Input genbank file." )

    # output file
    requiredNamed.add_argument( '--outfile', '-o', metavar = 'FILENAME', type = argparse.FileType( 'w', encoding = 'UTF-8' ),required = True, help = 'Name of output file (e.g filename.txt or filename.tsv).' )

    args = parser.parse_args()




# output file
out_file = args.outfile

# header line for output table
header = "\t".join(['Locus_tag', 'Type', 'Contig', 'Contig_length', 'Start', 'End', 'Length', 'Strand', 'Inference', 'Note', 'Codon_start', 'Translation_table', 'Product', 'Protein_ID', 'Translation', 'DNA_sequence', "\n"])

out_file.write(header)



# get all sequence records for the genbank file
recs = [ rec for rec in SeqIO.parse( args.infile, "genbank" ) ]

# iterate over each sequence record (contig in draft genome) and extract info,
# skipping features source, gene, operon
skip = [ 'source', 'gene', 'operon' ]
for rec in recs:

# NOTE that Main features, eg  accession No, etc. are NOT provided in Prokka annotations
    print( len(recs) )                # number of extracted contigs or sequence records
    organism = rec.description          # DEFINITION field
    if 'accessions' in rec.annotations.keys():
        acc = rec.annotations['accessions']
    else:
        acc = 'NA'

    if 'sequence_version' in rec.annotations.keys():
        ver = rec.annotations['sequence_version']
    else:
        ver = 'NA'

    acc_ver = '.'.join([acc, ver])      # VERSION field - most recent accession number

    print( organism )
    print( acc_ver )

    contig = rec.name # OR rec.id       # LOCUS field
    contig_len = len( rec.seq )         # contig length

    if rec.features:
        # get all the features in a given contig skipping source, gene, operon
        feats = rec.features

        for feat in feats:
            if feat.type not in skip:
                seq = feat.location.extract(rec).seq # sequence
                strand = str(feat.strand)
                if strand == '-1': # rev_com if opposite strand
                    seq = seq.reverse_complement()
                seq = str(seq)
                start = int( feat.location.nofuzzy_start ) # start, end coordinates
                end = int( feat.location.nofuzzy_end )
                length = len(seq)

                # feature type
                f_type = feat.type # feature type (e.g. CDS, rRNA, etc.)

                # not all features have all qualifiers so
                # if they're absent add 'NA'
                if 'locus_tag' in feat.qualifiers.keys():
                    locus = feat.qualifiers['locus_tag'][0] # locus_tag
                else:
                    locus = 'NA'

                if 'inference' in feat.qualifiers.keys():
                    inference = feat.qualifiers['inference'][0] # annotation tool
                else:
                    inference = 'NA'

                if 'note' in feat.qualifiers.keys():
                    note = feat.qualifiers['note'][0]
                else:
                    note = 'NA'

                if 'codon_start' in feat.qualifiers.keys():
                    codon_start = feat.qualifiers['codon_start'][0]
                else:
                    codon_start = 'NA'

                if 'transl_table' in feat.qualifiers.keys():
                    transl_table = feat.qualifiers['transl_table'][0] # translation table
                else:
                    transl_table = 'NA'

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


                # join features into tab-separated line
                line_out = "\t".join([ locus, f_type, contig, str(contig_len), str(start), str(end), str(length), strand, inference, note, str(codon_start), str(transl_table), product, protein_id, translation, seq, "\n" ])
                # print to file
                out_file.write(line_out)

out_file.close()

