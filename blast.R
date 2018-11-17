# http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# blastn -db db/db_dgri -query Drosophila_melanogaster.BDGP6.cdna.all.fa \
# -out results.out -outfmt 10
prefix <- 'gi/assignment_2/'
cols <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
             'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
blastn <- read.csv(paste0(prefix, 'dgri_blastn.txt'), col.names=cols)
tblastn <- read.csv(paste0(prefix, 'dgri_tblastn.txt'), col.names=cols)

# total blast alignments
print(paste('Total blastn alignments:', nrow(blastn)))
print(paste('Total tblastn alignments:', nrow(tblastn)))

# total unique genes found
blastn.uniq <- unique(blastn$qseqid)
tblastn.uniq <- unique(tblastn$qseqid)

length(blastn.uniq) # 13511 uniq qid
length(tblastn.uniq) # 29887 uniq qid

# bin data 
bins <- function(blast.data, thresh1, thresh2) {
  bin1 <- blast.data$evalue <= thresh1
  bin2 <- blast.data$evalue > thresh1 & blast.data$evalue <= thresh2
  bin3 <- blast.data$evalue >= thresh2
  x <- list(bin1, bin2, bin3)
  # lapply(x, function(bin) mean(bin))
  list(bin1=blast.data[bin1,], bin2=blast.data[bin2,], bin3=blast.data[bin3,])
}

thresh1 <- 1E-50
thresh2 <- 0.01

# filter with threshold
df <- blastn[blastn$evalue <= thresh1,]
test <- c('FBtr0089289', 'FBtr0089289', 'FBtr0100448', 'FBtr0089290', 'FBtr0089288', 'FBtr0100447', 'FBtr0329933', 'FBtr0336949')

matches <- vector()
for (query in unique(blastn$qseqid)) {
  # query <- 'FBtr0089289'
  hits <- df[df$qseqid == query,]
  uniq.scaff <- unique(hits$sseqid)
  matches <- c(matches, length(uniq.scaff))
}

hist(matches)

matches <- vector()
for (seq in unique(blastn$sseqid)) {
  hits <- df[df$sseqid == seq,]
  uniq.query <- unique(hits$qseqid)
  matches <- c(matches, length(uniq.query))
}

plot.ecdf(matches)


# find all unique cdna matches, assuming only match once.
# explore outlier cdna that match 10+ times

matches <- vector()
for (query in unique(blastn$qseqid)) {
  hits <- df[df$qseqid == query,]
  uniq.scaff <- unique(hits$sseqid)
  if (length(uniq.scaff) > 1) {
    matches <- c(matches, query)
  }
}
length(matches)
