rand_genomics_utils
======

*Utility scripts for random common tasks in genomics* 

Description
-----------

Repository of scripts without a common theme, for performing various tasks in comparative genomics. 


#### Script genbank_to_tsv.py

Python script to parse genbank record and convert to tab-delimited file.
Tested on bacterial genbank records, both complete and draft.
Feature qualifiers specific to eukaryotic records are not captured here but can 
easily be added.


```
./genbank_to_tsv.py -h
usage: genbank_to_tsv.py [-h] --infile GENBANK --outfile FILENAME

Convert genbank record into tab-separated table.

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  --infile GENBANK, -i GENBANK
                        Input genbank file.
  --outfile FILENAME, -o FILENAME
                        Name of output file (e.g filename.txt or
                        filename.tsv).
```


#### Script SNVs_from_alignmnt2.pl

Perl script to parse SNVs from multiple sequence alignment (MSA) of genes from e.g. clustalo or your favorite gene MSA algorithm.  This script expects the MSA to be FASTA formatted. The script parses nucleotide sequence alignments of genes to output tab-delimited SNVs and nucleotide positions for each sequence along with amino acid translations showing synonymous or non-synonymous changes.  It randomly takes the first sequence in the MSA as the reference and SNVs in all other sequences are called with respect to that reference sequence.  The script will take 


```
Usage:
SNVs_from_alignmnt2.pl MSA.fasta > MSA_SNVs.tsv
```
#### Script gbk2fasta_protein.py

Python script to extract gene amino acid translations from GENBANK file and format as FASTA.  Header contains gene locus_tag, protein ID, and functional characterization.


```
Usage:
gbk2fasta_protein.py filename.gbk > filename.faa
```
