
# parse genscan output

#scf_oneHit <- readLines(con="genscan/gs_scaffold_4820.fasta.out")

#scf_manyHits <- readLines(con="genscan/gs_scaffold_6488.fasta.out")

#scf_noHits <- readLines(con = "genscan/gs_scaffold_5604.fasta.out")

parseAndReturnList <- function(x){
  # parses the genscan output and returns a list with 2 elements,
  # first element is name of scaffold and info
  # second element is dataframe of found proteins
  
  #x <- scf_manyHits
  nameOfOutput <- x[[3]]
  
  # checks in no exons found
  if("NO EXONS/GENES PREDICTED IN SEQUENCE" %in% x){
    return(list(nameOfOutput, NA))
  }
  
  lastLineOfTable <- (match("Predicted peptide sequence(s):", x) - 2)
  
  #extracts only information to be used in table
  chrOfHits <- x[c(9, 12:lastLineOfTable)]
  
  chr_clean <- chrOfHits[chrOfHits != ""] 
  
  # c(1:5), c(7:10), c(12), c(14:19), c(21:26), c(28:31), c(33:34), c(36:37), c(39:42), c(44:47), c(49:53), c(55:59), c(61:66)
  
  # incides of text that is useful
  # all genscan tables have 66 char per table row
  indicesOfGenscanTable <- list(c(1:5), c(7:10), c(12), c(14:19), c(21:26), c(28:31), c(33:34), c(36:37), c(39:42), c(44:47), c(49:53), c(55:59), c(61:66))
  
  # preallocate for result
  listOfRows <- vector(mode = "list", length = length(chr_clean)+1)
  
  # go through every line of input and return only text from the indices that are relevant
  for(i in 1:length(chr_clean)){
    #i <- 1
    returnText <- lapply(indicesOfGenscanTable, function(indices, tableLine){
      #indices <- indicesOfGenscanTable[[1]]
      #tableLine <- find_40
      
      splitted <- strsplit(tableLine, "")[[1]]
      chars <- splitted[indices]
      chars_noEmpty <- chars[chars != " "]
      return(paste0(chars_noEmpty, collapse = ""))
      
    }, tableLine = chr_clean[i])
    listOfRows[[i]] <- returnText
  }
  
  # store in dataframe
  mtrx <- matrix(unlist(listOfRows), ncol=13, byrow=TRUE)
  df <- as.data.frame(mtrx, stringsAsFactors = F)
  # make colnames out of first row
  colnames(df) <- df[1,]
  df <- df[-1,]
  df2 <- as.data.frame(lapply(df, type.convert), stringsAsFactors = F)
  
  return(list(nameOfOutput, df2))
}

#t1 <- parseAndReturnList(scf_noHits)
#t2 <- parseAndReturnList(scf_oneHit)
#t3 <- parseAndReturnList(scf_manyHits)

# read every file in directory
# I had all genscan results in the directory genscan and they started with gs_
allFiles <- list.files("genscan/", pattern = "gs_*")

# for every file
for(i in 1:length(allFiles)){
  #i <- 1
  oneGenscanOutput <- readLines(con = paste0("genscan/", allFiles[i]))
  
  # parse the file
  parsed <- parseAndReturnList(oneGenscanOutput)
  
  # write csv file
  outputFile <- paste0("genscan/", allFiles[i], ".csv")
  
  write(parsed[[1]], file =  outputFile)
  
  if(is.na(parsed[[2]])){
    write(parsed[[2]], file = outputFile, append = T)
  }else{
    write.table(parsed[[2]], file =  outputFile, append = TRUE, 
                quote = FALSE, sep = ",", na = "NA", dec = ".",
                col.names = TRUE, row.names = F)
  }
  
}


