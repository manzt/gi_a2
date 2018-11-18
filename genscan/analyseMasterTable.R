# analyse genscan

library(ggplot2)

analyseMasterTable <- function(mt){
  #mt <- masterTable
  colsSplit <- strsplit(mt$Gn.Ex, ".", fixed = T)
  gn <- lapply(colsSplit, function(x) x[1])
  ex <- lapply(colsSplit, function(x) x[2])
  
  masterTableGE <- cbind(mt, gene = (unlist(gn)), exon = (unlist(ex)))
  
  masterTableGE$gene <- as.character(levels(masterTableGE$gene)[masterTableGE$gene])
  masterTableGE$exon <- as.character(levels(masterTableGE$exon)[masterTableGE$exon])
  
  #   add score label to score
  #   about score: a predicted 5' splice site with score > 100 is strong; 50-100 is moderate;
  #   0-50 is weak; and below 0 is poor (more than likely not a real donor site).
  #   range_ -20 : 2000
  scoresOrdered <- data.frame(matrix(data = c( -30, -1, 0, 49, 50, 99, 100, 2500),
                                     ncol = 2, byrow = T))
  colnames(scoresOrdered) <- c("min", "max")
  
  colWithScore <- sapply(masterTableGE$Tscr.., function(score){
    l <- NA
    if(!is.na(score)){
      l <- which(score >= scoresOrdered[,"min"] & score <= scoresOrdered[,"max"])
    }
    return(l)
  })
  scoreLabel <- unlist(as.character(colWithScore))
  scoreLabelP <- as.numeric(scoreLabel)
  
  masterTableGES <- cbind(masterTableGE, scoreLabelP)
  
  # number of scaffolds with no predicted genes
  
  scaffolds <- split(masterTableGES, f = masterTableGES$name)
  
  genesFound <- lapply(scaffolds, function(x){
    #x <- scaffolds[["scaffold_1009"]]
    if(is.na(x[1,"Gn.Ex"])){
      return(FALSE)
    }else{
      return(TRUE)
    }
  })
  
  scaffoldContainsGene <- sum(unlist(genesFound))
  scaffoldContainsNoGene <- sum(!unlist(genesFound))
  
  # number of genes per scaffold
  scaffoldsWithGenes <- scaffolds[unlist(genesFound)]
  
  numOfGenes <- lapply(scaffoldsWithGenes, function(x){
    #x <- scaffoldsWithGenes[[1]]
    numOfG <- max(as.numeric(x$gene))
  })
  
  # mean
  meanGenesPerScf <- mean(unlist(numOfGenes))
  # hist
  hist(unlist(numOfGenes))
  # total genes predicted
  sum(unlist(numOfGenes))
  
  # play with probability values
  pDist <- lapply(scaffoldsWithGenes, function(a){
    return(a[,"P...."])
  })
  
  probVals <- as.numeric(unlist(pDist, use.names = F))
  probVals <- probVals[!is.na(probVals)]
  df_probVals <- as.data.frame(cbind(probVals, idx = c(1:length(probVals))))
  
  plottyplot <- ggplot(df_probVals, aes(x = probVals)) + geom_histogram()
  
  
  # remove all predicted genes where one or more exons have a probability of less than 0.5
  numOfGenesWithGoodExonsProb <- lapply(scaffoldsWithGenes, function(x){
    #x <- scaffoldsWithGenes[[1]]
    scfByGene <- split(x, f = x$gene)
    geneGood <- lapply(scfByGene, function(y){
      #y <- scfByGene[[1]]
      ps <- y$P....[which(!is.na(y$P....))]
      if(any(ps < 0.5) | length(ps) == 0){
        return(FALSE)
      }else{
        return(TRUE)
      }
    })
    return(sum(unlist(geneGood)))
  })
  
  # sum(unlist(numOfGenesWithGoodExonsProb))
  # 18412

  returnGenesAboveBoundary <- function(x, boundary){
    #x <- scaffoldsWithGenes[[1]]
    scfByGene <- split(x, f = x$gene)
    
    geneGoodScore <- lapply(scfByGene, function(y){
      #y <- scfByGene[[3]]
      scores <- y$scoreLabel[which(!is.na(y$scoreLabel))]
      if(any(scores < boundary) | length(scores) == 0){
        return(FALSE)
      }else{
        return(TRUE)
      }
    })
    return(sum(unlist(geneGoodScore)))
  }
    
  numOfGenesWithGoodExonsScore <- lapply(scaffoldsWithGenes, returnGenesAboveBoundary,
                                         boundary = 2)
  sum(unlist(numOfGenesWithGoodExonsScore))
  
  
  
  #
  
  # number of exons per gene per scaffold
  
  # 
  
  
  
  
  
  
  
}







