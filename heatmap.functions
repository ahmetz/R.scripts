rpkm.heatmap <- function(dataframe, title){
  require("RColorBrewer")
  dmatrix <- data.matrix(dataframe[2:length(dataframe)])
  dmatrix <- apply(dmatrix, 2, function(x) log2(x+1))
  row.names(dmatrix) <- dataframe$gene
  hr <- hclust(as.dist(1-cor(t(dmatrix), method="pearson")), method="complete")
  hc <- hclust(as.dist(1-cor(dmatrix, method="spearman")), method="complete")
  y<-dmatrix
  nf <- layout
  heatmap(dmatrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale = "none", col= brewer.pal(9, "Blues"), main=title, margins = c(10,10))

  heatmap(matrix(rep(seq(min(as.vector(y)), max(as.vector(y)), length.out=10),2), 2, byrow=T, dimnames=list(1:2, round(seq(min(as.vector(y)), max(as.vector(y)), length.out=10),1))), col=brewer.pal(9, "Blues"), Colv=NA, Rowv=NA, labRow="", main=title)
}
getRPKMdata <- function(geneList){
  connection <- dbConnect(SQLite(), "~/sqlDatabase/geneExpression.db") #get connected to the SQL database
  temp2 <- geneList[,1] # convert gene names into a vector that is iterable
  result <- data.frame() # intialize the data frame to store the RPKM data
  notFound <- data.frame()
  #geneList <- as.vector(geneList) #make sure the list is a vector
  for(i in (1:length(temp2))){ #start the loop for reading each gene and getting RPKM values from the SQL database
    gene <- temp2[i]
    query <- paste('SELECT * FROM tissue_rpkm WHERE gene="', gene,'"', sep="") # form the query
    temp <- dbGetQuery(connection, query) #store results of the query in temp
    if(length(row.names(temp))!=0){
      temp.cast <- cast(temp, gene ~ tissue, value = "RPKM") #cast the data in a long format
      result <- rbind(result, temp.cast)
    }else{
    notFound <- c(notFound, gene)
    }}

  return(result)
  rpkm.heatmap(result, title)
  if(length(notFound)==0){
    paste("All genes have been found in the database")
  }else{
  paste("These genes do not have a corresponding gene symbol in the database", notFound)
  }
}

genesToHeatmap <- function(geneList){
  dataframe <- getRPKMdata(geneList)
  title <- deparse(substitute(geneList))
  rpkm.heatmap(dataframe, title)

  }

