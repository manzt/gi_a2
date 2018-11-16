# split me
#install.packages("seqinr")
library(seqinr)

scaffolds <- seqinr::read.fasta(file = "dmoj_all_chromosome.fasta")

lengthsOfScaf <- lapply(scaffolds, function(x){

  return(length(x))

})

vecLen <- unlist(lengthsOfScaf, use.names = F)

# check lengths of scaffolds
vecLen[1:20]>500000

# cutoff at ~7 million bp
# for scaffolds 1:5 (which are all longer than 7mil bp): separate with
# overlap and write fasta files

trash <- lapply(scaffolds[1:12], function(x){
  
  #x <- scaffolds[[1]]
  parts <- ceiling((length(x))/500000)
  
  startHere <- 1
  endHere <- 600000
  for(i in 1:parts){
    #cat(i, "\n")
    write.fasta((x[startHere:endHere]), names = attributes(x)$name, file.out = paste0("genscan/", attributes(x)$name, "_part_", i, ".fasta"))
    startHere <- startHere + 500000
    endHere <- endHere + 500000
    if(endHere > length(x)){
      endHere <- length(x)
    }
    
  }
  return(NA)
})

# creates 325 files


# write one fasta file per scaffold
trash <- lapply(scaffolds[13:6841], function(x){
  #x <- scaffolds[[11]]
  write.fasta(x, names = attributes(x)$name, file.out = paste0("genscan/", attributes(x)$name, ".fasta"))
  return(NA)
})

# creates  files
# in total: 7154 files

# now run genscan on all these files
