library(ggplot2)
library(dplyr)
library(UpSetR)

url.prefix <- 'https://raw.githubusercontent.com/manzt/gi_a2/master/all-transcripts-results/'
# url.prefix <- 'gi/assignment_2/'
species <- c('Dgrim', 'Dmoj', 'Dvir')
# species <- c('dgri')
cols <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
          'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

raw_data <- list()
for (s in species) {
  url <- paste0(url.prefix, s, 'Blast.txt')
  print(url)
  raw_data[[s]] <- read.csv(url, col.names=cols) %>% mutate(species=s)
}
blastn <- bind_rows(raw_data)

# total blast alignments for each species
blastn %>% 
  group_by(species) %>% 
  tally()

ggplot(blastn, aes(x=species)) + geom_bar()

# define e-value threshold and filter data
# thresh <- 1E-50
thresh <- 800
blastn <- blastn %>% filter(bitscore < thresh)

# total blast alignments for each species after thresh filter
blastn %>% 
  group_by(species) %>% 
  tally()

ggplot(blastn, aes(x=species)) + geom_bar()


# count transcripts found from Dmel 
found_genes <- blastn %>%
  group_by(species, sseqid) %>%
  summarise()

# number of Dmel genes in each
found_genes %>% tally()

# plot Intersection
by_species <- lapply(species, function(s) found_genes$sseqid[found_genes$species == s])
names(by_species) <- species
upset(fromList(by_species), order.by='freq')

