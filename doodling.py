# import LethrinidaeSpecies.csv
import csv
import numpy as np
from Bio import Entrez, SeqIO, AlignIO
from Bio.Nexus import Nexus
import time
import os
import sys
import General as g
import glob
from collections import Counter
import re
import pickle

Entrez.email = 'hannahiweller@gmail.com'

"""EVENTUAL WORKFLOW GOAL:
1) Input SPECIES LIST and GENE LIST (i.e. "COI", "16S", "Cytb", "rag1", etc)
2) Search for accession numbers of each species/gene combination
3) Store the accession numbers by gene
4) Download all Genbank files for stored accession numbers by gene
5) Scrub (quality control) sequences for < 0.03 proportion ambiguous characters and longer than minimum length for that gene
6) Convert scrubbed GBK files to FASTA for review in Mesquite (pick only 1 sequence/species for every gene)**
7) Review aligned Nexus files in Mesquite and save final version (aligned, 1 sequence/species for every gene) in SIMPLIFIED Nexus format
8) Edit taxon names down to just species names instead of huge Genbank tags
9) Align all genes into supermatrix

** For Sparids, also then tack this onto the bottom of Lethrinids gene matrices and save as simplified nexus."""

# input SPECIES LIST and GENE LIST
projectName = 'Lutjanidae' # project name
speciesList = g.flatten(list(zip(*csv.reader(open('CSV/Lutjanidae.csv', 'rb'))))) # opens 1D species list (from rfishbase)
del speciesList[0] # remove column name ('Species' is first element of list when you flatten it)

geneList = ["COI", "16S", "CytB", "RAG1", "Rhodopsin"]




# path to CSV with genbank accession numbers
# LethrinidPath = '/Users/hannah/Dropbox/Westneat_Lab/GenBank Scraping/LethrinidaeSpecies.csv'
# projectName = 'Lethrinidae'



# reader = csv.reader(open('Sparids.csv', 'rb'))
# byColumn = list(zip(*reader))
# genes = ["COI", "16S", "Cytb", "rag1", "Rhodopsin"]
# genes = ["Cytb", "rag1", "Rhodopsin"]

# for gene in genes:
#     record = []
#     saveName = 'Sparidae' + gene + 'AccNumbers'
#     for species in byColumn[0]:
#         searchTerm = species + "[organism] AND " + gene + "[title]"
#         handle = Entrez.esearch(db="nucleotide", term=searchTerm)
#         temp = Entrez.read(handle)
#         if temp["Count"] != '0':
#             record.append(temp)
#     pickle.dump(record, open(saveName, 'w'))

# pickle.dump(record, open('./SparidAccNumbers', 'w'))
#
# geneAccNumbers = glob.glob('./Sparidae*AccNumbers')

# for gene in geneAccNumbers:
#     SparidAccNumbers = pickle.load(open(gene, 'rb'))
#     saveName = gene[0:-1-9]+'Raw.gbk'
#     gbkNumbers = []
#     for entry in SparidAccNumbers:
#         gbkNumbers.append(entry["IdList"])
#
#     gbkNumbers = g.flatten(gbkNumbers)
#     g.gbFetch(gbkNumbers, saveName)

# rawGBK = glob.glob('./*Raw.gbk')
# simplified = glob.glob("NEX/Lethrinidae*Simple.nex")
#
# combos = glob.glob("NEX/COMBINED*Simple.nex")
#
# nexi = [(fname, Nexus.Nexus(fname)) for fname in combos]
#
# combined = Nexus.combine(nexi)
# combined.write_nexus_data(filename=open('COMBINED_LethrinidSparid_All.nex', 'w'))

# for f in simplified:
#     input_handle = AlignIO.read(open(f), "nexus")
#     saveName = f[0:(-1-12)] + "_Renamed.nex"
#     output_handle = open(saveName, "w")
#
#     for f in input_handle:
#         newName = " ".join(f.id.split(" ")[1:3])
#         f.id = newName
#         print f.id
#
#     AlignIO.write(input_handle, output_handle, "nexus")

# g.scrubGBK(projectName='Sparidae', gbkFiles = rawGBK, geneNames = genes)
# fileNames = glob.glob('./*Scrubbed.gbk')
# g.fastaSeq(fileNames)
# scrubGBK = glob.glob('./*Scrubbed.gbk')
