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
import pickle

Entrez.email = 'hannahiweller@gmail.com'

# path to CSV with genbank accession numbers
LethrinidPath = '/Users/hannah/Dropbox/Westneat_Lab/GenBank Scraping/LethrinidaeSpecies.csv'
projectName = 'Lethrinidae'


# fetches all the sequences from the provided accession numbers and saves one gbk file per gene/column
# g.fetchSeq(LethrinidPath, projectName)
#
# # checks for length and ambiguous characters
# g.scrubSeq(LethrinidPath, projectName)
#
# # converts to fasta format
# gbkFiles = glob.glob('./*Scrubbed.gbk')
# g.fastaSeq(gbkFiles)

# pull out one sequence for each species - longest?

# how to pull down every sequence for a given species
# Entrez.email = 'hannahiweller@gmail.com'
# handle = Entrez.esearch(db="nucleotide", term="Monotaxis grandoculis[PORG]", rettype="gb")
# record = Entrez.read(handle)

# reader = csv.reader(open(LethrinidPath, 'rb')) # read in csv
# byColumn = list(zip(*reader)) # parse by columns, first element is definition
# Lethrinids = byColumn[0][1:-1]
#
# # g.getSaveGenBankIDs(Lethrinids, 'LethrinidAccessionNumbers')
#
# with open('LethrinidAccessionNumbers', 'rb') as f:
#     record = pickle.load(f)
#
# for species in record:
#     gbID = species["IdList"]
#
#     g.gbFetch(gbID, 'output.gbk')
#     time.sleep(1)

reader = csv.reader(open('Sparids.csv', 'rb'))
byColumn = list(zip(*reader))

record = []
# for species in byColumn[0]:
#     for gene in genes:
#         searchTerm = species + "[organism] AND " + gene + "[gene]"
#         handle = Entrez.esearch(db="nucleotide", term=searchTerm)
#         record.append(Entrez.read(handle))

genes = ["COI", "16S", "Cytb", "rag1", "Rhodopsin"]

for species in byColumn[0][1:-1]:
    searchTerm = species + "[organism] AND rag1[gene]"
    handle = Entrez.esearch(db="nucleotide", term=searchTerm)
    record.append(Entrez.read(handle))
# accNumbers = []
#
# for species in Lethrinids:
#     accNumbers.append(g.speciesSeq(species))
#
# with open('temp', 'wb') as f:
#     pickle.dump(accNumbers, f)

# with open('temp', 'rb') as f:
#     test = pickle.load(f)
# pickle.dump(accNumbers, 'temp')

# itemList = pickle.load('temp')




# handle = Entrez.einfo(db="nucleotide")
# record = Entrez.read(handle)




# for field in record["DbInfo"]["FieldList"]:
#     print("%(Name)s, %(FullName)s, %(Description)s" % field)

# PORG, Primary Organism, Scientific and common names of primary organism, and all higher levels of taxonomy
