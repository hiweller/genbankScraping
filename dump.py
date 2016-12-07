# A TRASH PLACE FOR DOODLE CODE WHEN I THINK I'M DONE WITH IT
# BUT NOT SO CONVINCED I'M DONE THAT I'M WILLING TO DELETE IT


# input SPECIES LIST and GENE LIST
# projectName = 'Lutjanidae' # project name
# speciesList = g.flatten(list(zip(*csv.reader(open('CSV/Sparids.csv', 'rb'))))) # opens 1D species list (from rfishbase)
# del speciesList[0] # remove column name ('Species' is first element of list when you flatten it)
#
# geneList = ["COI", "16S", "CytB", "RAG1", "Rhodopsin"]
#
# # g.gbkFamilyGeneSearch('Nemipteridae', geneList, 'Accession Numbers/Nemipteridae')
#
# geneAccNumbers = glob.glob('Accession Numbers/Nemipteridae*')
#
# for gene in geneAccNumbers:
#     AccNumberSet = pickle.load(open(gene, 'rb'))
#     saveName = 'GBK'+gene[17:-1-9]+'Raw.gbk'
#     gbkNumbers = []
#     for entry in AccNumberSet:
#         gbkNumbers.append(entry["IdList"])
#
#     gbkNumbers = g.flatten(gbkNumbers)
#     g.gbFetch(gbkNumbers, saveName)
# for gene in geneAccNumbers:
#     AccNumberSet = pickle.load(open(gene, 'rb'))
#     saveName = gene[0:-1-9]+'Raw.gbk'
#     gbkNumbers = []
#     for entry in AccNumberSet:
#         gbkNumbers.append(entry["IdList"])
#
#     gbkNumbers = g.flatten(gbkNumbers)
#     g.gbFetch(gbkNumbers, saveName)


# g.gbkSpeciesGeneSearch(speciesList, geneList, 'Sparidae/Species')

# byGene = glob.glob('Accession Numbers/Gene*')
#
# bg = []
#
# for i in byGene:
#     x = 0
#     AccNumberSet = pickle.load(open(i, 'rb'))
#     for j in AccNumberSet:
#         x = x + int(j["Count"])
#     bg.append(x)
#
# byTitle = glob.glob('Accession Numbers/Title*')
#
# bt = []
#
# for i in byTitle:
#     x = 0
#     AccNumberSet = pickle.load(open(i, 'rb'))
#     for j in AccNumberSet:
#         x = x + int(j["Count"])
#     bt.append(x)
#
# bySpecies = glob.glob('Accession Numbers/Species*')
#
# bs = []
#
# for i in bySpecies:
#     x = 0
#     AccNumberSet = pickle.load(open(i, 'rb'))
#     for j in AccNumberSet:
#         x = x + int(j["Count"])
#     bs.append(x)
# # search for accession numbers of each species/gene combination
# g.gbkSpeciesGeneSearch(speciesList, geneList, 'Accession Numbers/Lutjanidae')

# dowload GBK files for stored accession numbers
# geneAccNumbers = glob.glob('Accession Numbers/Lutjanidae*AccNumbers')
#
# for gene in geneAccNumbers:
#     AccNumberSet = pickle.load(open(gene, 'rb'))
#     saveName = gene[0:-1-9]+'Raw.gbk'
#     gbkNumbers = []
#     for entry in AccNumberSet:
#         gbkNumbers.append(entry["IdList"])
#
#     gbkNumbers = g.flatten(gbkNumbers)
#     g.gbFetch(gbkNumbers, saveName)

# for gene in geneAccNumbers:
#     SparidAccNumbers = pickle.load(open(gene, 'rb'))
#     saveName = gene[0:-1-9]+'Raw.gbk'
#     gbkNumbers = []
#     for entry in SparidAccNumbers:
#         gbkNumbers.append(entry["IdList"])
#
#     gbkNumbers = g.flatten(gbkNumbers)
#     g.gbFetch(gbkNumbers, saveName)





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



# with open('temp', 'rb') as f:
#     record = pickle.load(f)
#
# for species in record:
#     print species["Term"]
#     print species["IdList"]






# LethrinidGenes = glob.glob("NEX/*Renamed.nex")
#
# nexi = [(fname, Nexus.Nexus(fname)) for fname in LethrinidGenes]
#
# combined = Nexus.combine(nexi)
# combined.write_nexus_data(filename=open('CombinedTest.nex', 'w'))






# simplified = glob.glob("NEX/*Simple.nex")

# f = "NEX/LethrinidaeRhodopsin_oneSimple.nex"
# input_handle = AlignIO.read(open(f), "nexus")
# saveName = f[0:(-1-12)] + "Renamed.nex"
# output_handle = open(saveName, "w")
#
# for f in input_handle:
#     newName = " ".join(f.id.split(" ")[1:3])
#     f.id = newName
#     print f.id
#
# AlignIO.write(input_handle, output_handle, "nexus")

# for f in simplified:
    # input_handle = AlignIO.read(open(f), "nexus")
    # saveName = f[0:(-1-12)] + "Renamed.nex"
    # output_handle = open(saveName, "w")
    #
    # for f in input_handle:
    #     newName = " ".join(f.id.split(" ")[1:3])
    #     f.id = newName
    #     print f.id
    #
    # AlignIO.write(input_handle, output_handle, "nexus")


# test = 'NEX/LethrinidaeCytB_one.nex'
# alignment = AlignIO.read(open("NEX/LethrinidaeRAG1_oneSimple.nex"), "nexus")
# for record in alignment:
#     newName =  " ".join(record.id.split(" ")[1:3])
#     record.id = newName
#     print record.id
#
# output_handle = open("NEX/test.nex", 'w')
# AlignIO.write(alignment, output_handle, "nexus")
# handle = open(test)
# for f in AlignIO.read(handle, 'nexus'):
#     print f




















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




#
# for species in byColumn[0][1:-1]:
#     for gene in genes:
#         searchTerm = species + "[organism] AND " + gene + "[gene]"
#         handle = Entrez.esearch(db="nucleotide", term=searchTerm)
#         record.append(Entrez.read(handle))
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
