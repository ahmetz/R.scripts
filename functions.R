
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
  
  
paste.data <- function (header = FALSE) {
    read.table(pipe("pbpaste"), header = header)
}

sparkLinePlot <- function(df, plot.file) {
    
    highest <- subset(df, outcomes == max(outcomes))
    lowest <- subset(df, outcomes == min(outcomes))
    
    p <- ggplot(df, aes(x=date, y=outcomes)) +
        geom_line() +
        opts(panel.border = theme_rect(linetype = 0),
             panel.background = theme_rect(colour = "white"),
             panel.grid.major = theme_blank(),
             panel.grid.minor = theme_blank(),
             axis.text.x = theme_blank(),
             axis.text.y = theme_blank(),
             axis.ticks = theme_blank()) +
                 ylab("") +
                 geom_point(data = lowest, size = 3, colour = alpha("red", 0.5)) +
                 geom_point(data = highest, size = 3, colour = alpha("blue", 0.5)) +
                 scale_y_continuous(formatter = comma) +
                 scale_x_date(name = "", major="months", minor="weeks", format="%b-%d")
    
    p
    
}

my_theme <- function(base_size = 18, base_family="Gill Sans Light", dkpanel=FALSE, legendPosition = "bottom", xAngle = 0) {
  thme <- theme(
  	  line = element_line(colour = "#838281", size = 0.5,
                          linetype = 1, lineend = "round"),
      rect = element_rect(fill = "white",
                          colour = NA, size = 0.5, linetype = 1),
      text = element_text(family = base_family, face = "plain",
                          colour = "black", size = base_size, hjust = 0.5, vjust = 0.5,
                          angle = 0, lineheight = 1),
      axis.text = element_text(size = rel(1), face = "italic"),
      axis.text.x = element_text(vjust = 1, angle = xAngle),
      axis.text.y = element_text(hjust = 1),
      axis.line = element_line(colour = "#d8d6d5"),
      axis.line.y = element_line(),
      axis.line.x = element_line(),
      axis.ticks = element_line(),
      axis.title = element_text(size = rel(1)),
      axis.title.x = element_text(),
      axis.title.y = element_text(angle = 90),
 
      axis.ticks.length = unit(base_size * 0.2 , "points"),
      axis.ticks.margin = unit(base_size * 0.5, "points"),
      
      legend.background = element_rect(color = "grey", linetype=1, size = rel(0.5)),
      legend.margin = unit(base_size * 2.5, "points"),
      legend.key = element_rect(),
      legend.key.size = unit(1.5, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = rel(0.8)),
      legend.text.align = NULL,
      legend.title = element_text(size = rel(1),  hjust = 1),
      legend.title.align = NULL,
      legend.position = legendPosition,
      legend.direction = NULL,
      legend.justification = "left",
      legend.box = "horizontal",
      
      panel.background = element_rect(linetype=0),
      panel.border = element_rect(fill = NA, color = NA),
      panel.grid.major = element_line(colour = "#d8d6d5", size=rel(1.25), linetype = 3),
      panel.grid.minor = element_line(colour = "#edecec", size = rel(1), linetype = 3),
      panel.margin = unit(c(1, 1, 1, 1), "lines"),
      
      strip.background = element_rect(fill = "white", colour = "#838281", linetype=1, size=1),
      strip.text = element_text(size = rel(0.8)),
      strip.text.x = element_text(),
      strip.text.y = element_text(),
      
      plot.background = element_rect(fill = "white", colour="white"),
      plot.title = element_text(face = "bold", size = rel(1.5), hjust=0.5),
      plot.margin = unit(c(6, 20, 9, 10) * 2, "points"),
      complete = TRUE)

  
  thme
}

heatmap2 <- function(matrix){
  matrix <- as.matrix(matrix)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
 
hr <- hclust(as.dist(1-cor(t(matrix), method="pearson")), method="average")
  hc <- hclust(as.dist(1-cor(matrix, method="spearman")), method="average")
  heatmap.2(matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myPalette(100), scale="row", density.info="none", trace="none", keysize = 0.7)
}

TCGA.survival.p <- function(geneName, stratify = "mRNA"){
  #This function returns the p value for survival analysis for a fiven gene in TCGA database for ovarian cancer
  # if the input is a gene name then it will call the TCGA.gene function to obtain the data frame (TCGA object) which has the expression and survival data
  if(class(geneName) =="character") {
    gene <- TCGA.gene(geneName)
  } else {
    gene <- geneName
  }
  # gene <- data.frame(OS = gene$OS, OS.status = gene$OS.status, mRNA = gene$mRNA)
  gene$status <- gene$OS.status == "DECEASED"
  gene$mRNAlevels <- gene[stratify] > colMeans(gene[stratify])
  survival <- Surv(gene$OS, gene$status)
  lin.fit <- survfit(survival ~ gene$mRNAlevels)
  sdf <- survdiff(survival ~ gene$mRNAlevels)
  p <- 1-pchisq(sdf$chisq, length(sdf$n) -1)
  p
}

TCGA.survival.fit <- function(gene){
  #this function returns the survival fit object for a given TCGA object. Object should have survival and expression data
  # patients will be stratified based on expression mean values
  gene <- data.frame(OS = gene$OS, OS.status = gene$OS.status, mRNA = gene$mRNA)
  gene$status <- gene$OS.status == "DECEASED"
  gene$mRNAlevels <- gene$mRNA > mean(gene$mRNA)
  survival <- Surv(gene$OS, gene$status)
  lin.fit <- survfit(survival ~ gene$mRNAlevels)
  sdf <- survdiff(survival ~ gene$mRNAlevels)
  p <- 1-pchisq(sdf$chisq, length(sdf$n) -1)
  lin.fit
}

TCGA.RRPA <- function(geneName){
  # This function returns the RPPA data for a given gene from the TCGA ovarian cancer project
  # If the input is a gene, it will create the TCGA object by calling upon TCGA.gene function, otherwise it expects a TCGA gene object as input
  if(class(geneName) =="character") {
    gene <- TCGA.gene(geneName)
  } else {
    gene <- geneName
  }
  gene["RPPA"]
}

TCGA.gene <- function(geneName){
  # Main function to create a TCGA gene object. Input is a gene name. 
  # It does not yet check to make sure the gene name is in TCGA database yet
  # returns a dataframe, what I'd like to call a TCGA object
  library("cgdsr")
  library(ggplot2)
  library(xtable)
  library("gplots")
  library("reshape")
  library(survival)
  options(digits = 2)
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  geneexp1 <- getProfileData(mycgds, geneName, c("ov_tcga_pub_mrna_median", "ov_tcga_pub_rae", "ov_tcga_pub_methylation"), "ov_tcga_pub_exp1")
  geneexp2 <- getProfileData(mycgds, geneName, c("ov_tcga_pub_mrna_median", "ov_tcga_pub_rae", "ov_tcga_pub_methylation"), "ov_tcga_pub_exp2")
  geneexp3 <- getProfileData(mycgds, geneName, c("ov_tcga_pub_mrna_median", "ov_tcga_pub_rae", "ov_tcga_pub_methylation"), "ov_tcga_pub_exp3")
  geneexp4 <- getProfileData(mycgds, geneName, c("ov_tcga_pub_mrna_median", "ov_tcga_pub_rae", "ov_tcga_pub_methylation"), "ov_tcga_pub_exp4")
  geneexp1$patient <- row.names(geneexp1)
  geneexp2$patient <- row.names(geneexp2)
  geneexp3$patient <- row.names(geneexp3)
  geneexp4$patient <- row.names(geneexp4)
  geneexp1$cluster <- "Fallopian"
  geneexp2$cluster <- "Immunoreactive"
  geneexp3$cluster <- "Mesenchymal"
  geneexp4$cluster <- "Proliferative"
  row.names(geneexp1) <- NULL
  row.names(geneexp2) <- NULL
  row.names(geneexp3) <- NULL
  row.names(geneexp4) <- NULL
  geneAll <- data.frame()
  geneAll <- rbind(geneexp1, geneexp2)
  geneAll <- rbind(geneAll, geneexp3)
  geneAll <- rbind(geneAll, geneexp4)
  clinical <- getClinicalData(mycgds, "ov_tcga_pub_all")
  clinical$patient <- row.names(clinical)
  row.names(clinical) <- NULL
  geneAll <- join(geneAll, clinical, type="inner")
  gene <- data.frame(patient = geneAll$patient, cluster = geneAll$cluster, mRNA = geneAll$ov_tcga_pub_mrna_median, RAE = geneAll$ov_tcga_pub_rae, meth = geneAll$ov_tcga_pub_methylation, OS = geneAll$overall_survival_months, OS.status = geneAll$overall_survival_status)
  gene
}

TCGA <- function(geneName, expGraph = TRUE, singleSurvGraph = TRUE, clusterSurvGraph = TRUE){
  gene <- TCGA.gene(geneName)
  print(paste("------- Gene Name:", geneName, "------------"))
  if(sum(is.na(gene$mRNA == 489))) {
    print(paste("There is no mRNA expression available for", geneName))
    stop("There is no mRNA expression available for this gene")
  }
  gene$status <- gene$OS.status == "DECEASED"
  gene$mRNAlevels <- gene$mRNA > mean(gene$mRNA)
  
  a <- ggplot(data = gene, aes(x = factor(cluster), y = mRNA)) + geom_boxplot() + theme_bw(18) + xlab("") + ylab("Normalized expression") + labs(title=paste("Expression of" , geneName, "in TCGA database"))
  par(mfrow=c(3, 2), mar = c(3, 3, 3, 3))
  if(expGraph == TRUE) with(gene, plotmeans(mRNA ~ cluster, xlab="Expression clusters based on Nature paper", ylab=paste(geneName, "Expression"), main=paste(geneName, "expression in 4 different \ncluster groups of TCGA patients")))
  m <- ifelse(mean(gene$mRNA) > gene$mRNA, 0, 1)
  percent <- sum(m)/length(m)*100
  print(sprintf("The percentage of patients with expression levels above the mean  is: %.1f", percent))
  p <- TCGA.survival.p(gene)
  print(sprintf("The survival p value is %.2f", p))
  if(singleSurvGraph == TRUE) {
    lin.fit <- TCGA.survival.fit(gene)
    plot(lin.fit, lty=1:2, mark.time=F, bg="white")
    legend(10, .2, c(paste(geneName, "High Expression"), paste(geneName, "Low Expression")), lty=1:2, cex = 0.8)
    title(main=paste("Survival Probablity of \nTCGA Ovarian Patients \n based on ", geneName, " expression (p=", format(p, digits = 2),  ")", sep = ""), xlab="Months", ylab="Probability")
  }
  
  b <- split(gene, gene$cluster)
  c <- sapply(b, TCGA.survival.p)
  d <- sapply(b, TCGA.survival.fit)
  if(clusterSurvGraph == TRUE) {
  #  par(mfrow=c(2, 2), mar = c(2, 2, 2, 2))
    plot(d[1]$Fallopian, lty=1:2, mark.time=F, bg="white", main=paste("Fallopian p=", format(c[1] , digits=2)))
    legend(10, .2, c("High", "Low"), lty=1:2, cex = 0.8)
    plot(d[2]$Immunoreactive, lty=1:2, mark.time=F, bg="white", main = paste("Immunoreactive p=", format(c[2] , digits=2)))
    legend(10, .2, c("High", "Low"), lty=1:2, cex = 0.8)
    plot(d[3]$Mesenchymal, lty=1:2, mark.time=F, bg="white",  main = paste("Mesenchymal p=", format(c[3] , digits=2)))
    legend(10, .2, c("High", "Low"), lty=1:2, cex = 0.8)
    plot(d[4]$Proliferative, lty=1:2, mark.time=F, bg="white", main = paste("Proliferative p=", format(c[4] , digits=2)))
    legend(10, .2, c("High", "Low"), lty=1:2, cex = 0.8)
  }
  data.frame(geneName, expressionPercent = percent, Overall.p.value = p, Fallopian.p.value = c[1], Immunoreactive.p.value = c[2], Mesenchymal.p.value = c[3], Proliferative.p.value = c[4])
}

geneSymbolsToEntrez <- function(geneSymbols) {
  library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
getBM(attributes="entrezgene", filters = "hgnc_symbol", values=geneSymbols, mart = ensembl)

}

TCGA.expression <- function(genes){#get the expression data from TCGA database
  library(cgdsr)
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  expression <- data.frame(rep(1, 563)) #Create a dataframe to put the date in
  for (i in 1:length(genes)) { #chr.list is a list of genes
    x <- getProfileData(mycgds, genes[i], "ov_tcga_pub_mrna_median", "ov_tcga_pub_all")
    if(length(table(x)) != 0 ) {  #Check to make sure the gene is represented and expression values are obtained
      expression <- cbind(expression, x)
    }
  } 
  df <- na.omit(expression[2:length(genes)]) #Remove missing values
  df
}

survivalMultipleGenes <- function(genes, title){
  df <- TCGA.expression(genes)
  df$patient <- rownames(df) # add a column with patient IDs
  df.melt <- melt(df, id = "patient") # melt the data
  df.melt <- ddply(df.melt, .(variable), transform, mean = mean(value)) # create a new column, with mean value for each gene in the group of patients
  df.melt <- ddply(df.melt, .(patient), transform, PIR = ifelse(value > mean, 1, -1)) #For each gene, if the expression is higher than the mean, then the PIR score is 1, else it is -1
  df.PIR <- ddply(df.melt, .(patient), colwise(sum, "PIR")) # Sum the PIR values
  df.PIR <- ddply(df.PIR, .(patient), transform, survival = ifelse(PIR > median(df.PIR$PIR), 1, 0)) #stratify patients based on the PIR value
  clinical <- getClinicalData(mycgds, "ov_tcga_pub_all") #Get clinical data
  clinical$patient <- rownames(clinical)
  df.PIR <- join(df.PIR, clinical, type = "inner")
  df.PIR$status <- df.PIR$overall_survival_status=="DECEASED"
  df.surv <- Surv(df.PIR$overall_survival_months, df.PIR$status)
  df.diff <- survdiff(df.surv ~ df.PIR$survival)
  df.fit <- survfit(df.surv ~ df.PIR$survival)
  p <- 1-pchisq(df.diff$chisq, length(df.diff$n) -1)
  plot(df.fit, lty=1:2, mark.time=F)
  legend(10, .2, c("High", "low"), lty = 1:2)
  title(main=paste("Survival Probablity of TCGA Ovarian Patients \nbased on ", title,"  expression (p = ", round(p, digits = 2), ")", sep=""), xlab="Months", ylab="Probability")
}

GSEAtoGraph <- function(GSEAfile1, GSEAfile2, sampleName1 = NULL, sampleName2 = NULL){
  a <- read.delim(GSEAfile1, nrows = 10)
  a <- a[, c(1, 4, 6)]
  a$NAME <- gsub("_", " ", a$NAME)
  plot1 <- ggplot(a, aes(x = NES, y = NAME, size = SIZE)) + geom_point(color = "red") + my_theme(base_size = 14, base_family = "Helvetica") + labs(title = paste("Gene set enrichment in genes upregulated in\n", sampleName1, " vs ", sampleName2), x = "Normalized Enrichment Score", y = " ") + scale_size_area(max_size=10, name = "Number of genes belonging to a gene set")
  b <- read.delim(GSEAfile2, nrows = 10)
  b <- b[, c(1, 4, 6)]
  b$NAME <- gsub("_", " ",b$NAME)
  plot2 <- ggplot(b, aes(x = NES, y = NAME, size = SIZE)) + geom_point(color = "red") + my_theme(base_size = 14, base_family = "Helvetica") + labs(title = paste("Gene set enrichment in genes upregulated in\n", sampleName2, " vs ", sampleName1), x = "Normalized Enrichment Score", y = " ") + scale_size_area(max_size=10, name = "Number of genes belonging to a gene set")
  multiplot(plot1, plot2, cols = 2)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##########################################
## Intersect and Venn Diagram Functions ##
##########################################
## Author: Thomas Girke
## Last update: March 24, 2012
## Utilities: 
## (1) Venn Intersects
##     Computation of Venn intersects among 2-20 or more sample sets using the typical
##     'only in' intersect logic of Venn comparisons, such as: objects present only in 
##     set A, objects present only in the intersect of A & B, etc. Due to this restrictive 
##     intersect logic, the combined Venn sets contain no duplicates.  
## (2) Regular Intersects
##     Computation of regular intersects among 2-20 or more sample sets using the
##     following intersect logic: objects present in the intersect of A & B, objects present 
##     in the intersect of A & B & C, etc. The approach results usually in many duplications 
##     of objects among the intersect sets.
## (3) Graphical Utilities
##     - Venn diagrams of 2-5 sample sets. 
##     - Bar plots for the results of Venn intersect and all intersect approaches derived 
##       from many samples sets. 
##
## Detailed instructions for using the functions of this script are available on this page:
##     http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn
##
## Revision history:
##     March 24, 2012: fixed substring problem in plotVenn function

#######################################
## Define Generic Intersect Function ##
#######################################
## Computation of (1) Venn Intersects and (2) Regular Intersects
overLapper <- function(setlist=setlist, complexity=1:length(setlist), sep="-", cleanup=FALSE, keepdups=FALSE, type) {
  ## Clean up of sample sets to minimize formatting issues 
  if(cleanup==TRUE) {
    ## Set all characters to upper case 
    setlist <- sapply(setlist, function(x) gsub("([A-Z])", "\\U\\1", x, perl=T, ignore.case=T))
    ## Remove leading and trailing spaces
    setlist <- sapply(setlist, function(x) gsub("^ {1,}| {1,}$", "", x, perl=T, ignore.case=T))
  }
  
  ## Append object counter to retain duplicates 
  if(keepdups==TRUE) {
    dupCount <- function(setlist=setlist) {
      count <- table(setlist)
      paste(rep(names(count), count), unlist(sapply(count, function(x) seq(1, x))), sep=".")
    }
    mynames <- names(setlist)
    setlist <- lapply(setlist, function(x) dupCount(x)) # lapply necessary for numeric data!
    names(setlist) <- mynames
  }	
  
  ## Create intersect matrix (removes duplicates!)
  setunion <- sort(unique(unlist(setlist)))
  setmatrix <- sapply(names(setlist), function(x) setunion %in% unique(setlist[[x]])) 
  rownames(setmatrix) <- setunion
  storage.mode(setmatrix) <- "numeric"
  
  ## Create all possible sample combinations within requested complexity levels
  labels <- names(setlist)
  allcombl <- lapply(complexity, function(x) combn(labels, m=x, simplify=FALSE))
  allcombl <- unlist(allcombl, recursive=FALSE)
  complevels <- sapply(allcombl, length)
  
  ## Return intersect list for generated sample combinations 
  if(type=="intersects") {
    OLlist <- sapply(seq(along=allcombl), function(x) setunion[rowSums(setmatrix[, rep(allcombl[[x]], 2)]) == 2 * length(allcombl[[x]])])
    names(OLlist) <- sapply(allcombl, paste, collapse=sep)
    return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Intersect_List=OLlist))
  }	
  
  ## Return Venn intersect list for generated sample combinations 
  if(type=="vennsets") {
    vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
      mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
      mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
      cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
      cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
      return(setunion[cond1 & cond2])
    }
    vennOLlist <- sapply(seq(along=allcombl), function(x) vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x))
    names(vennOLlist) <- sapply(allcombl, paste, collapse=sep)
    return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Venn_List=vennOLlist))
  }
}

###########################################
## Define Venn Diagram Plotting Function ##
###########################################
vennPlot <- function(counts=counts, mymain="Venn Diagram", mysub="default", setlabels="default", yoffset=seq(0,10,by=0.75), ccol=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), colmode=1, lcol=c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00"), lines=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"), mylwd=3, diacol=1, type="ellipse", ccex=1.0, lcex=1.0, sepsplit="_", ...) {
  library(RColorBrewer)
  ## Enforce list structure to support multiple venn sets 
  if(is.list(counts)==FALSE) {
    counts <- list(counts)
  }
  
  ## Check for supported number of Venn counts: 3, 7, 15 and 31
  if(!length(counts[[1]]) %in%  c(3,7,15,31)) stop("Only the counts from 2-5 way venn comparisons are supported.")
  
  ## Function to return for a set label the index of matches in the name field of a counts object
  grepLabel <- function(label, x=names(counts[[1]])) {
    x <- strsplit(x, sepsplit)
    as.numeric(which(sapply(x, function(y) any(y==label))))
  }
  
  ## 2-way Venn diagram
  if(length(counts[[1]])==3) {
    ## Define subtitle
    if(mysub=="default") {
      n <- names(counts[[1]])[1:2]
      if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
      } else {
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
      }
      mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), sep="")
    } else { 
      mysub <- mysub 
    }
    
    ## Plot venn shapes
    symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
    
    ## Add counts
    for(i in seq(along=counts)) {
      olDF <- data.frame(x=c(3.1, 7.0, 5.0), 
                         y=c(6.0, 6.0, 6.0), 
                         counts=counts[[i]])
      if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
      if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) } # For coloring several numbers per intersect differently. ccol can needs to be list to color each field differently..
    }
    
    ## Add sample labels
    if(length(setlabels)==1 & setlabels[1]=="default") { 
      setlabels <- names(counts[[1]][1:2])
    } else {
      setlabels <- setlabels
    }
    text(c(1.0, 9.0), c(8.8, 8.8), labels=setlabels, col=lcol, cex=lcex, ...)	
  }
  
  ## 3-way Venn diagram
  if(length(counts[[1]])==7) { 
    ## Define subtitle
    if(mysub=="default") {
      n <- names(counts[[1]])[1:3]
      if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
      } else {
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
      }
      mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), sep="")
    } else { 
      mysub <- mysub
    }
    
    ## Plot venn shapes
    symbols(x=c(4, 6, 5), y=c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)
    
    ## Add counts
    for(i in seq(along=counts)) {
      olDF <- data.frame(x=c(3.0, 7.0, 5.0, 5.0, 3.8, 6.3, 5.0), 
                         y=c(6.5, 6.5, 3.0, 7.0, 4.6, 4.6, 5.3), 
                         counts=counts[[i]])
      if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
      if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
      
    }
    
    ## Add sample labels
    if(length(setlabels)==1 & setlabels[1]=="default") { 
      setlabels <- names(counts[[1]][1:3])
    } else {
      setlabels <- setlabels
    }
    text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels=setlabels, col=lcol, cex=lcex, ...)	
  }
  
  ## 4-way Venn diagram with ellipses
  if(length(counts[[1]])==15 & type=="ellipse") {
    ## Define subtitle
    if(mysub=="default") {
      n <- names(counts[[1]])[1:4]
      if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
      } else {
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
      }
      mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
    } else { 
      mysub <- mysub
    }
    
    ## Plot ellipse
    plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
      angles <- (0:segments) * 2 * pi/segments  
      rotate <- rotate*pi/180
      ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
      ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
      ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
      plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
    }
    ## Plot ellipse as 4-way venn diagram
    ellipseVenn <- function(...) {
      split.screen(c(1,1))
      plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
      screen(1, new=FALSE)
      plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
      screen(1, new=FALSE)
      plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
      screen(1, new=FALSE)
      plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
      ## Add counts
      for(i in seq(along=counts)) {
        olDF <- data.frame(x=c(1.5, 3.5, 6.5, 8.5, 2.9, 3.1, 5.0, 5.0, 6.9, 7.1, 3.6, 5.8, 4.2, 6.4, 5.0), 
                           y=c(4.8, 7.2, 7.2, 4.8, 5.9, 2.2, 0.7, 6.0, 2.2, 5.9, 4.0, 1.4, 1.4, 4.0, 2.8), 
                           counts=counts[[i]])
        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
        if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
      }
      ## Add sample labels
      if(length(setlabels)==1 & setlabels[1]=="default") { 
        setlabels <- names(counts[[1]][1:4])
      } else {
        setlabels <- setlabels
      }
      text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels=setlabels, col=lcol, cex=lcex, ...)
      close.screen(all=TRUE) 
    }
    ellipseVenn(...)
  } 
  
  ## 4-way Venn diagram with circles (pseudo-venn diagram that misses two overlap sectors) 
  if(length(counts[[1]])==15 & type=="circle") {
    ## Define subtitle
    if(mysub=="default") {
      n <- names(counts[[1]])[1:4]
      if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
      } else {
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
      }
      mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
    } else { 
      mysub <- mysub
    }
    
    ## Plot venn shapes
    symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)
    
    ## Add counts
    for(i in seq(along=counts)) {
      olDF <- data.frame(x=c(3.0, 6.5, 3.0, 6.5, 4.8, 3.0, 4.8, 4.8, 6.5, 4.8, 3.9, 5.7, 3.9, 5.7, 4.8), 
                         y=c(7.2, 7.2, 3.2, 3.2, 7.2, 5.2, 0.4, 0.4, 5.2, 3.2, 6.3, 6.3, 4.2, 4.2, 5.2), 
                         counts=counts[[i]])
      if(colmode==1) { text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol, cex=ccex, ...) } # rows 14-15 of olDF are printed in next step
      if(colmode==2) { text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol[[i]], cex=ccex[i], ...) }
      text(c(4.8), c(0.8) + yoffset[i], paste("Only in ", names(counts[[1]][1]), " & ", names(counts[[1]][4]), ": ", olDF$counts[7], "; Only in ", names(counts[[1]][2]), " & ", names(counts[[1]][3]), ": ", olDF$counts[8], sep=""), col=diacol, cex=ccex, ...)
    }
    
    ## Add sample labels
    if(length(setlabels)==1 & setlabels[1]=="default") { 
      setlabels <- names(counts[[1]][1:4])
    } else {
      setlabels <- setlabels
    }
    text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels=setlabels, col=lcol, cex=lcex, ...)
  } 
  
  ## 5-way Venn diagram
  if(length(counts[[1]])==31) {
    ## Define subtitle
    if(mysub=="default") {
      n <- names(counts[[1]])[1:5]
      if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
      } else {
        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
      }
      mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), paste("; S5 =", sample_counts[5]), sep="")
    } else { 
      mysub <- mysub
    }
    
    ## Plot ellipse
    plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
      angles <- (0:segments) * 2 * pi/segments  
      rotate <- rotate*pi/180
      ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
      ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
      ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
      plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
    }
    ## Plot ellipse as 5-way venn diagram
    ellipseVenn <- function(...) {
      split.screen(c(1,1))
      screen(1, new=FALSE)
      plotellipse(center=c(4.83,6.2), radius=c(1.43,4.11), rotate=0, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
      screen(1, new=FALSE)
      plotellipse(center=c(6.25,5.4), radius=c(1.7,3.6), rotate=66, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
      screen(1, new=FALSE)
      plotellipse(center=c(6.1,3.5), radius=c(1.55,3.9), rotate=150, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
      screen(1, new=FALSE)
      plotellipse(center=c(4.48,3.15), radius=c(1.55,3.92), rotate=210, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
      screen(1, new=FALSE)
      plotellipse(center=c(3.7,4.8), radius=c(1.7,3.6), rotate=293.5, segments=360, xlab="", ylab="", col=lines[5], axes=FALSE, lwd=mylwd, ...)
      
      ## Add counts
      for(i in seq(along=counts)) {
        olDF <- data.frame(x=c(4.85, 8.0, 7.1, 3.5, 2.0, 5.90, 4.4, 4.60, 3.60, 7.1, 6.5, 3.2, 5.4, 6.65, 3.40, 5.00, 6.02, 3.60, 5.20, 4.03, 4.20, 6.45, 6.8, 3.39, 6.03, 5.74, 4.15, 3.95, 5.2, 6.40, 5.1), 
                           y=c(8.30, 6.2, 1.9, 1.6, 5.4, 6.85, 6.6, 2.45, 6.40, 4.3, 6.0, 4.6, 2.1, 3.40, 3.25, 6.43, 6.38, 5.10, 2.49, 6.25, 3.08, 5.30, 4.0, 3.80, 3.20, 5.95, 5.75, 3.75, 3.0, 4.50, 4.6),
                           counts=counts[[i]]) 
        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
        if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
      }
      ## Add sample labels
      if(length(setlabels)==1 & setlabels[1]=="default") { 
        setlabels <- names(counts[[1]][1:5])
      } else {
        setlabels <- setlabels
      }
      text(c(5.7, 7.9, 8.5, 4.2, 0.8), c(9.9, 7.9, 1.9, 0.0, 7.3), adj=c(0, 0.5), labels=setlabels, col=lcol, cex=lcex, ...)
      close.screen(all=TRUE) 
    }
    ellipseVenn(...)
  } 
}

##############################
## Define Bar Plot Function ##
##############################
## Plots the counts of Venn/regular intersects generated by the overLapper function
olBarplot <- function(OLlist=OLlist, mycol="default", margins=c(6, 10, 3, 2), mincount=0, mysub="default", ...) {
  ## Generate counts and allow lower limit 
  counts <- sapply(OLlist[[4]], length)
  mylogical <- counts >= mincount
  counts <- counts[mylogical]
  
  ## Color bars by default by complexity levels 
  if(mycol=="default") {
    mycol <- OLlist$Complexity_Levels
    mycol <- mycol[mylogical] 
  } else {
    mycol <- mycol	
  }
  
  ## Define subtitle
  if(mysub=="default") {
    mysub <- paste("Min Count:", mincount)
  } else {
    mysub <- mysub
  }
  
  ## Generate bar plot with defined margins
  par(mar=margins) # Define margins to allow long labels
  barplot(counts, col=mycol, sub=mysub, ...)
  par(mar=c(5, 4, 4, 2) + 0.1) # Set margins back to default
}


