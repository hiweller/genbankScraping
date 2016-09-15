# PRESETS/RANGE RESTRICTIONS FOR COMMONLY ENCOUNTERED BARCODE GENES
import csv
import numpy as np
from Bio import Entrez, SeqIO
import time
import os
import sys
import re
import glob
import pickle

lengthRestrictions = {'12S':500, '16S':525, 'RAG1':1200, 'RAG2':1200, 'RHODOPSIN':500, 'CYTB':700, 'COI':500}

AmbiguousChars = ['R', 'Y', 'W', 'S', 'M', 'K', 'H', 'B', 'D', 'V', 'X', 'N']
PropAmThresh = 0.03
GBMax = 100

# function for splitting into chunks
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def fetchSeq(csvPath, projectName=''):
    reader = csv.reader(open(csvPath, 'rb')) # read in csv
    byColumn = list(zip(*reader)) # parse by columns, first element is definition

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
        gbFetch(gene, filename)

def scrubSeq(csvPath, projectName=''):
    """Sorting FASTA sequences with following conditions:
    1) Remove anything with > 0.03 ambiguous character proportion
    2) Remove anything shorter than recommended minimum gene length for that gene."""
    gbkFiles = glob.glob('./*.gbk')
    reader = csv.reader(open(csvPath, 'rb'))
    geneNames = next(reader)[1:len(next(reader))]

    for gene in gbkFiles:

        KeepPile = []
        MaybePile = []
        ThrowawayPile = []

        searchTerms = re.compile(r"(?=("+'|'.join(geneNames)+r"))", re.IGNORECASE)
        geneName = re.findall(searchTerms, gene)[0]

        geneList = lengthRestrictions.keys()
        LR = None

        saveName = projectName + geneName + 'Scrubbed.gbk'

        if geneName.upper() in map(str.upper, lengthRestrictions.keys()):
            LR = lengthRestrictions[geneName.upper()]

        for f in SeqIO.parse(gene, 'genbank'):
            seq = f.seq
            species = f.annotations['organism']
            geneName = f.description

            PropAmbig = sum(seq.count(x) for x in AmbiguousChars)/len(seq)

            if PropAmbig < PropAmThresh:
                MaybePile.append(f)
            else:
                ThrowawayPile.append(f)

        for f in MaybePile:
            if LR is not None and len(f.seq) < LR:
                ThrowawayPile.append(f)
            else:
                KeepPile.append(f)

        if len(KeepPile) is not 0:
            output_handle = open(saveName, 'w')
            count = SeqIO.write(KeepPile, output_handle, 'genbank')

def gbFetch(IDs, filename): # pass list of GB id's
    if len(IDs) > GBMax:
        for chunk in chunks(IDs[1:len(IDs)], 100):
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
        genbank_handle = Entrez.efetch(db='nucleotide', id=IDs[1:len(IDs)], rettype='gb', retmode='text')
        if os.path.isfile(filename):
            out_handle = open(filename, 'a')
        else:
            out_handle = open(filename, 'w')
        out_handle.write(genbank_handle.read())
        out_handle.close()
        genbank_handle.close()

def getSaveGenBankIDs(species, filename):
    """Takes list of species, searches GenBank for vouchers associated with each and saves them."""
    output = []

    for s in species:
        output.append(speciesSeq(s))

    with open(filename, 'wb') as f:
        pickle.dump(output, f)

def speciesSeq(species):
    searchTerm = species + "[PORG]"
    Entrez.email = 'hannahiweller@gmail.com'
    handle = Entrez.esearch(db="nucleotide", term=searchTerm, rettype="gb")
    record = Entrez.read(handle)
    return record

def fastaSeq(fileNames):
    for gbk in fileNames:
        savename = gbk[0:len(gbk)-4]+'.fasta'
        input_handle = open(gbk, 'rU')
        output_handle = open(savename, 'w')

        sequences = SeqIO.parse(input_handle, 'genbank')
        count = SeqIO.write(sequences, output_handle, 'fasta')

        output_handle.close()
        input_handle.close()
