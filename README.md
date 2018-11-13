# gi_a2
Genome Informatics, Assignment 2: Repo contains outputs for blast alignments as well as R scripts.

### BLAST output files
The headers for the (t)blastn file output files are as follows:
`'qid', 'sid', 'pctIdnt', 'length', 'qstart', 'qend', 'sstart', 'send', 'eval', 'score'`


### File naming

#### blastn

| Query        | Database  | File name  |
| ------------- |:-------------:| :-----:|
| cDNA Dmel  | Dvir all-chrom  | `dvir_blastn.txt` |
| - - | Dgri all-chrom   |   `dgri_blastn.txt` |
|- - | Dmaj all-chrom      |   `dmaj_blastn.txt` |

#### tblastn

| Query  | Database | File name |
| ------------- |:-------------:| :-----:|
| pep Dmel  | Dvir all-chrom  | `dvir_tblastn.txt` |
| - -  | Dgri all-chrom  | `dgri_tblastn.txt` |
| - -  | Dmaj all-chrom  | `dmaj_tblastn.txt`  |
