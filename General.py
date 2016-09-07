# PRESETS/RANGE RESTRICTIONS FOR COMMONLY ENCOUNTERED BARCODE GENES
# 12S, 16S swapped to S12, S16 to avoid initializing variable w/ numbers
import csv
import numpy as np
from Bio import Entrez
import time
import os
import sys

S12 = 500
S16 = 525
RAG1 = 1200
RAG2 = 1200
Rhodopsin = 500
CytB = 700
COI = 500
AmbiguousChars = ['R', 'Y', 'W', 'S', 'M', 'K', 'H', 'B', 'D', 'V', 'X', 'N']
PropAmThresh = 0.03

# function for splitting into chunks
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def fetchSeq(csvPath, projectName=''):
    reader = csv.reader(open(csvPath, 'rb')) # read in csv
    byColumn = list(zip(*reader)) # parse by columns, first element is definition

    GBMax = 100

    GBCodes = []
    for el in byColumn:
        el = [i.split(',') for i in el]
        el = sum(el, [])
        el = filter(None, el)
        GBCodes.append(el)

    for gene in GBCodes[1:len(GBCodes)]:

        # get gene name and build filename
        filename = projectName + gene[0] + '.gbk'

        # fetching FASTA sequences
        # if more than 100, space it out!
        if len(gene) > GBMax:
            for chunk in chunks(gene[1:len(gene)], 100):
                genbank_handle = Entrez.efetch(db='nucleotide', id=chunk, rettype='gb', retmode='text')
                if os.path.isfile(filename):
                    out_handle = open(filename, 'a')
                    out_handle.write(genbank_handle.read())
                else:
                    out_handle = open(filename, 'w')
                    out_handle.write(genbank_handle.read())
                out_handle.close()
                genbank_handle.close()
                time.sleep(3)
        else:
            genbank_handle = Entrez.efetch(db='nucleotide', id=gene[1:len(gene)], rettype='gb', retmode='text')
            out_handle = open(filename, 'w')
            out_handle.write(genbank_handle.read())
            out_handle.close()
            genbank_handle.close()
