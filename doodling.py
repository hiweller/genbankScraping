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
import re

Entrez.email = 'hannahiweller@gmail.com'

# path to CSV with genbank accession numbers
LethrinidPath = '/Users/hannah/Dropbox/Westneat_Lab/GenBank Scraping/LethrinidaeSpecies.csv'
projectName = 'Lethrinidae'


# fetches all the sequences from the provided accession numbers and saves one gbk file per gene/column
g.fetchSeq(LethrinidPath, projectName)

# checks for length and ambiguous characters
g.scrubSeq(LethrinidPath, projectName)

# converts to fasta format
gbkFiles = glob.glob('./*Scrubbed.gbk')
g.fastaSeq(gbkFiles)







# fastaFiles = glob.glob('./*.fasta')
# gbkFiles = glob.glob('./*.gbk')


# reader = csv.reader(open(LethrinidPath, 'rb'))
# geneNames = next(reader)[1:len(next(reader))]
#
# for gene in gbkFiles:
#
#     KeepPile = []
#     MaybePile = []
#     ThrowawayPile = []
#
#     searchTerms = re.compile(r"(?=("+'|'.join(geneNames)+r"))", re.IGNORECASE)
#     geneName = re.findall(searchTerms, gene)[0]
#
#     geneList = g.lengthRestrictions.keys()
#     LR = None
#
#     saveName = projectName + geneName + 'Scrubbed.gbk'
#
#     if geneName.upper() in map(str.upper, g.lengthRestrictions.keys()):
#         LR = g.lengthRestrictions[geneName.upper()]
#
#     for f in SeqIO.parse(gene, 'genbank'):
#         seq = f.seq
#         species = f.annotations['organism']
#         geneName = f.description
#
#         PropAmbig = sum(seq.count(x) for x in g.AmbiguousChars)/len(seq)
#
#         if PropAmbig < g.PropAmThresh:
#             MaybePile.append(f)
#         else:
#             ThrowawayPile.append(f)
#
#     for f in MaybePile:
#         if LR is not None and len(f.seq) < LR:
#             ThrowawayPile.append(f)
#         else:
#             KeepPile.append(f)
#
#     if len(KeepPile) is not 0:
#         output_handle = open(saveName, 'w')
#         count = SeqIO.write(KeepPile, output_handle, 'genbank')
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
