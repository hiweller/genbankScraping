# import LethrinidaeSpecies.csv
import csv
import numpy as np
from Bio import Entrez, SeqIO, AlignIO
from Bio.Nexus import Nexus
from Bio import Seq
import time
import os
import sys
import General as g
import glob
from collections import Counter
import re
import pickle
from Bio.Align.Applications import MuscleCommandline

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

"""Include:
1) Lethrinidae
2) Sparidae
3) Nemipteridae
4) Sillaginidae
5) Lobotidae
6) Scatophagidae
7) Siganidae
8) Caproidae
9) Moronidae
10) Lutjanidae
11) Serranidae
12) Lophiiformes"""

familyList = ["Nemipteridae", "Sillaginidae", "Lobotidae", "Scatophagidae", "Siganidae", "Caproidae", "Moronidae", "Lutjanidae", "Serranidae"]

geneList = ["COI", "16S", "CytB", "RAG1", "Rhodopsin"]

# grabbing accession numbers for genes in listed families:
# for fam in familyList:

    # # we'll store all collected data files (pickle, FASTA, GBK, Nexus) in a folder named for the family we're searching
    # saveFolder = "./"+fam
    # if not os.path.exists(saveFolder):
    #     os.makedirs(saveFolder)
    #
    # # for redundancy, name each saved file with family name + gene name
    # saveName = saveFolder + "/" + fam
    #
    # # saves accession numbers by gene for all species in family in pickled python files called familyNameGeneNameAccNumbers in respective family folders
    # g.gbkFamilyGeneSearch(fam, geneList, saveName)

# take accession numbers and download associated sequences:
# results in 'Raw' and 'Scrubbed' (QC) GBK files
# scrubbed are also converted to FASTA files
# for fam in familyList:
#
#     # get directory for accession numbers + read them in
#     folderName = fam+'/*AccNumbers'
#     geneAccNumbers = glob.glob(folderName)
#
#     for gene in geneAccNumbers:
#         # load each gene accession number set
#         AccNumberSet = pickle.load(open(gene, 'rb'))
#
#         # save in same folder with family name + gene name again
#         # 'Raw' because there will also be a QC version in a second
#         saveName = fam + '/' + fam + gene[17:-1-9]+'Raw.gbk'
#         gbkNumbers = []
#
#         # download GBK files for each one + save in *Raw.gbk files
#         for entry in AccNumberSet:
#             gbkNumbers.append(entry["IdList"])
#
#         gbkNumbers = g.flatten(gbkNumbers)
#         g.gbFetch(gbkNumbers, saveName)
#
#     # scrub for ambiguous characters or too-short sequences
#     g.scrubSeq(filePath = fam, geneList = geneList, projectName = fam+'/'+fam)
#
#     fileNames = glob.glob(fam + "/" + fam + "*Scrubbed.gbk")
#     g.fastaSeq(fileNames)
#
#     # convert to FASTA files for viewing in Mesquite
#
# # do muscle alignment on each fasta file
# for fam in familyList:
#     fastaFiles = glob.glob(fam+'/*.fasta')
#     muscle_cmd = "/Users/hannah/Dropbox/Westneat_Lab/Phylo Resources/muscle3.8.31_i86darwin64"
#     for fas in fastaFiles:
#         in_file = fas
#         out_file = fas[0:-6]+'.afa'
#         muscle_cline = MuscleCommandline(muscle_cmd, input = in_file, out = out_file)
#         muscle_cline()

# alphabetizing by taxon...
# for fam in familyList:
#     AFAs = glob.glob('Temp output/'+fam+'/*.afa')
#
#     for a in AFAs:
#         sequence = AlignIO.read(open(a, 'rU'), 'fasta')
#         sortedSequence = [f for f in sorted(sequence, key=lambda x : " ".join(x.description.split(" ")[1:3]))]
#
#         saveName = a[0:(-1-3)]+'Sorted.afa'
#         output_handle =  open(saveName, 'w')
#         SeqIO.write(sortedSequence, output_handle, 'fasta')

## INTERMISSION ##
# open in Mesquite + select sequences you actually want

# combine the lethrinid+sparid matrices with appropriate matrices from all other families
combined = glob.glob('NEX/COMBINED_*Simple.nex')
addFish = glob.glob('Temp output/Simplified Nexus/*Simple.nex')

# for f in combined:
#
#     input_handle = AlignIO.read(open(f), "nexus")
#     records = list(input_handle)
#
#     # maxLen = len(input_handle[0])# max(len(record.seq) for record in records)
#
#     searchTerms = re.compile(r"(?=("+'|'.join(geneList)+r"))", re.IGNORECASE)
#     geneName = re.findall(searchTerms, f)[0]
#
#     for fam in addFish:
#         if geneName in fam:
#             add_handle = AlignIO.read(open(fam), "nexus")
#             records2 = list(add_handle)
#             maxLen = len(input_handle[0])
#             maxLen2 = len(add_handle[0])#max(len(record.seq) for record in records2)
#
#             if maxLen < maxLen2:
#                 for sequence in input_handle:
#                     newSeq = str(sequence.seq).ljust(maxLen2, '-')
#                     sequence.seq = Seq.Seq(newSeq)
#             elif maxLen2 <= maxLen:
#                 for sequence in add_handle:
#                     newSeq = str(sequence.seq).ljust(maxLen, '-')
#                     sequence.seq = Seq.Seq(newSeq)
#
#             for line in add_handle:
#                 input_handle.append(line)
#
#
#     saveName = 'Test/'+geneName+'.nex'
#     output_handle = open(saveName, 'w')
#     AlignIO.write(input_handle, output_handle, "nexus")

lethrinids = glob.glob('NEX/COMBINED*Renamed_Simple.nex')
outgroups = glob.glob('Outgroups/*Simple.nex')

nexi = [(j, Nexus.Nexus(j)) for j in lethrinids]
lethrinidsCombined = Nexus.combine(nexi)
lethrinidsCombined.write_nexus_data(filename=open('OUT/Lethrinids_combined_all.nex', 'w'))


nexi2 = [(j, Nexus.Nexus(j)) for j in outgroups]
outgroupsCombined = Nexus.combine(nexi2)
outgroupsCombined.write_nexus_data(filename=open('OUT/Outgroups_combined_all.nex', 'w'))


# combined = glob.glob('Test/*Simple.nex')
# nexi = [(j, Nexus.Nexus(j)) for j in combined]
# combined = Nexus.combine(nexi)
# combined.write_nexus_data(filename=open('Test/COMBINED_All.nex', 'w'))
            # combos = [f, fam]
            # nexi = [(j, Nexus.Nexus(j)) for j in combos]
            # combined = Nexus.combine(nexi)



# limits just to species name
# for f in addFish:
#     input_handle = AlignIO.read(open(f), "nexus")
#
#     for i in input_handle:
#         newName = " ".join(i.id.split(" ")[1:3])
#         i.id = newName
#         # i.description = newName
#         print i.id
#
#     # sortedSequence = [j for j in sorted(input_handle, key=lambda x : " ".join(x.description.split(" ")[1:3]))]
#
#     saveName = 'Temp output/Simplified Nexus/'+f.split("/")[-1]
#     output_handle = open(saveName, 'w')
#     SeqIO.write(input_handle, output_handle, "nexus")

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
