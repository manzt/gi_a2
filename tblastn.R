library(ggplot2)
library(dplyr)
library(UpSetR)

paths <- c('gi/assignment_2/dgri_tblastn.txt',
           '/mhome/maths/t/ldk26/private/GenomeInformatics/D.mojavensis/result_tblastn2.out',
           '/mhome/maths/s/as2841/GIA2/dvri_tblastn.txt')

species <- c('dgri', 'dmoj', 'dvir')

dgri.cols <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
          'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

dmoj.cols <- c('qseqid', 'sseqid', 'pident', 'length', 'qstart', 
               'sstart', 'send', 'evalue', 'bitscore')

dvir.cols <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

cols <- list(dgri.cols, dmoj.cols, dvir.cols)

raw_data <- list()
for (i in 1:3) {
  print(paths[i])
  raw_data[[species[i]]] <- read.csv(paths[i], col.names=cols[[i]]) %>% mutate(species=species[i])
  head(raw_data[[species[i]]])
}
tblastn <- bind_rows(raw_data)

# count tblastn hits
tblastn %>%
  group_by(species) %>%
  count()

# set threshold and filter
thresh <- 1E-50
filtered <- tblastn %>% filter(evalue < thresh)

# count filtered tblast hits
filtered %>%
  group_by(species) %>%
  tally()

ggplot(tblastn, aes(x=species)) + geom_bar()

# count transcripts found from Dmel 
found_genes <- tblastn %>%
  group_by(species, qseqid) %>%
  summarise()

# number of Dmel genes in each
found_genes %>% tally()

# plot Intersection
by_species <- lapply(species, function(s) found_genes$qseqid[found_genes$species == s])
names(by_species) <- species
upset(fromList(by_species), order.by='freq')


