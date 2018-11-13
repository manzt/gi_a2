# http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# blastn -db db/db_dgri -query Drosophila_melanogaster.BDGP6.cdna.all.fa \
# -out results.out -outfmt "10 qseqid sseqid pident length qstart qend sstart send eval bitscore"
filename <- 'dgri_blastn.txt'
HEADERS <- c('qid', 'sid', 'pctIdnt', 'length', 'qstart', 'qend', 'sstart', 'send', 'eval', 'score')
blastn.results <- read.csv(filename, col.names=HEADERS)

