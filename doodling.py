# import LethrinidaeSpecies.csv
import csv
import numpy as np
from Bio import Entrez, SeqIO
import time
import os
import sys
import General as g
import glob
from collections import Counter

Entrez.email = 'hannahiweller@gmail.com'

# path to CSV with genbank accession numbers
LethrinidPath = '/Users/hannah/Dropbox/Westneat_Lab/GenBank Scraping/LethrinidaeSpecies.csv'
projectName = 'Lethrinidae'

# fetches all the sequences from the provided accession numbers and saves one fasta file per gene/column
# g.fetchSeq(LethrinidPath, projectName)
fastaFiles = glob.glob('./*.fasta')
gbkFiles = glob.glob('./*.gbk')

"""Sorting FASTA sequences with following conditions:
1) Remove anything with > 0.03 ambiguous character proportion
2) Remove anything shorter than recommended minimum gene length for that gene."""


for gene in gbkFiles:

    KeepPile = []
    MaybePile = []
    ThrowawayPile = []

    for f in SeqIO.parse(gene, 'genbank'):
        seq = f.seq
        PropAmbig = sum(seq.count(x) for x in g.AmbiguousChars)/len(seq)

        if PropAmbig < g.PropAmThresh:
            MaybePile.append(f)
        else:
            ThrowawayPile.append(f)

        


#
# for gene in fastaFiles:
#     fasta_sequences = SeqIO.parse(open(gene), 'fasta')
#
#     for f in fasta_sequences:
#         seq = f.seq
#
#
# for fasta in fastaFiles:
#     fasta_sequences = SeqIO.parse(open(fasta), 'fasta')
#     for f in fasta_sequences:
#         print f
