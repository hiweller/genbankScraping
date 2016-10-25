setwd('/Users/hannah/Dropbox/Westneat_Lab/GenBank Scraping/CSV/')
library(rfishbase)
fishFam <- 'Lutjanidae'

speciesList <- species_list(Family = fishFam)

write.csv(x = data.frame(Species=speciesList), row.names=FALSE, file = paste(fishFam, '.csv', sep=""))
