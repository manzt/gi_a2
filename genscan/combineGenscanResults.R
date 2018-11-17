# combine all genscan outputs into 1 single matrix

allFilesCsv <- list.files("genscan/", pattern = "*.csv")

masterTable <- data.frame()

for(i in 1:length(allFilesCsv)){
  #i <- 5
  #inputFile <- "genscan/gs_scaffold_6473_part_1.fasta.out.csv"
  
  inputFile <- paste0("genscan/", allFilesCsv[i])
  firstLine <- readLines(con = inputFile, n = 1)
  tableOfExons <- read.table(file = inputFile, skip = 1, header = TRUE,
                             sep = ",")
  
  scfIdTrimBegin <- substring(inputFile, 12, nchar(inputFile))
  scfTrimEnd <- substring(scfIdTrimBegin, 1, nchar(scfIdTrimBegin) - 14)
  if(dim(tableOfExons)[1] == 0){
    
    asDf <- data.frame(matrix(NA,ncol = 13))
    colnames(asDf) <- c("Gn.Ex","Type","S",".Begin","...End",".Len","Fr","Ph",
                                      "I.Ac", "Do.T", "CodRg", "P....", "Tscr..")
    dataAll <- cbind(name = scfTrimEnd, description = firstLine, asDf)  
    
  }else{
    dataAll <- cbind(name = scfTrimEnd, description = firstLine, tableOfExons)
  }
  
  masterTable <- rbind(masterTable, dataAll)
  
}

#mastertable before: 181694
#sum of vecLen 193 826 310

