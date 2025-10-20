library(ggplot2)
require(ggalluvial)
library(randomcoloR)

# Create a dataframe for tracking clones
cdr3_dataframe.fx <- function(datapath, chain, filelist, totalinframe){
  
  if (!(totalinframe %in% c("total", "inframe"))) {
    stop("Error: unknown argument ", totalinframe, ". Please provide either total (for all clonotypes) or inframe (for in-frame clonotypes only)")
  }
  
  # Ensure only one chain is included
  filelist <- filelist[grepl(chain, filelist)]
  
  #Compile a big file with patient's mixcr files
  i <- 1
  for (f in filelist){
    mixcrfle <- read.table(paste(datapath, f, sep = ""), header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, na.strings = c("", "NA"))
    if(i == 1){
      compldfle <- mixcrfle
      compldfle <- cbind(cloneno = row.names(compldfle), filename = f, compldfle)
      i <- i + 1
    }
    else{
      compldfle1 <- mixcrfle
      compldfle1 <- cbind(cloneno = row.names(compldfle1), filename = f, compldfle1)
      compldfle <- rbind(compldfle, compldfle1)
      rm(compldfle1)
    }
  }
  myfiles <- unique(as.character(compldfle$filename))
  message("my files:")
  print(myfiles)
  
  message("Total recovered clonotypes:")
  print(length(compldfle$aaSeqCDR3))
  
  message("Total out-of-frame clonotypes:")
  print(length(compldfle$aaSeqCDR3[grepl("_", compldfle$aaSeqCDR3)]))
  message("Total clonotypes with stop codon:")
  print(length(compldfle$aaSeqCDR3[grepl("[*]", compldfle$aaSeqCDR3) &
                                     !grepl("_", compldfle$aaSeqCDR3)]))
  
  #make samplename column
  compldfle$filename <- as.character(compldfle$filename)
  compldfle$samplename <- gsub(".*.CLONES_","", compldfle$filename)
  compldfle$samplename <- gsub(chain,"", compldfle$samplename)
  
  # remove out-of-frame clonotypes and those with stop codon
  compldfle_clean <- compldfle[!grepl("_", compldfle$aaSeqCDR3) &
                                 !grepl("[*]", compldfle$aaSeqCDR3),]
  #Recalculate cloneFraction for each file
  compldfle_clean$cloneFraction <- NA
  for(f in myfiles){
    compldfle_clean$cloneFraction[compldfle_clean$filename == f] <- compldfle_clean$cloneCount[compldfle_clean$filename == f]/sum(compldfle_clean$cloneCount[compldfle_clean$filename == f])
  }
  
  message("Total productive clonotypes:")
  print(length(compldfle_clean$aaSeqCDR3))
  
  if(totalinframe == "inframe"){
    message("Output contains in_frame clonotypes only")
    return(compldfle_clean)}
  if(totalinframe == "total"){
    message("Output contains all clonotypes")
    return(compldfle)}
}

# Track a specific clone
track_Aclone.fx <- function(compldfle, plotpath, countfrac, clnefrc, clnseq){
  
  message("list of samples to track clones: ")
  mysamples <- unique(compldfle$samplename)
  print(mysamples)
  
  # Subset df
  CDR3_fraction <- compldfle[, c("samplename","nSeqCDR3","cloneFraction", "cloneCount")]
  
  # Subset to include only clonotypes with more than specified clonal fraction
  CDR3_fraction <- CDR3_fraction[CDR3_fraction$cloneFraction > clnefrc,]
  
  #Assign colors to the specific clone
  myclone <- clnseq
  notrecurring <- CDR3_fraction$nSeqCDR3[!CDR3_fraction$nSeqCDR3 %in% myclone]
  
  cloneColor <- distinctColorPalette(length(myclone))
  myColors <- c(cloneColor, rep("white",length(notrecurring)))
  names(myColors) <- c(myclone, notrecurring)
  
  # Generate a row for each sample that doesnot have jurkat clonotype
  ## This ensures alluvia are colored
  
  tmp <- CDR3_fraction[CDR3_fraction$nSeqCDR3 == myclone,]
  nonexisting <- mysamples[!mysamples %in% tmp$samplename]
  if(length(nonexisting) > 0){
    newentries <- data.frame("samplename" = nonexisting, "nSeqCDR3" = myclone,
                             "cloneFraction" = 0, "cloneCount" = 0)
    CDR3_fraction <- rbind(CDR3_fraction, newentries)
  }
  
  p <-  ggplot(CDR3_fraction, aes(x = samplename,
                                  y = eval(as.name(countfrac)),
                                  fill = nSeqCDR3,
                                  stratum = nSeqCDR3,
                                  alluvium = nSeqCDR3,
                                  label = nSeqCDR3))
  
  myp <- p + geom_alluvium(decreasing = FALSE) +
    geom_stratum(decreasing = FALSE, stat = "alluvium") +
    scale_fill_manual(breaks = names(myColors[myColors != "white"]),
                      values = myColors) +
    theme(axis.title.y = element_text(size = 50),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.position = "none",
          plot.margin = unit(c(0.2,0,0,0),"cm")) +
    labs(y = countfrac)
  
  pdf(paste0(plotpath, "trackAclone_", mysamples[1], countfrac, ".pdf"),
      width = 30,
      height = 20,
      useDingbats = FALSE,
      onefile = FALSE)
  print(myp)
  dev.off()
  
}





# For adaptive data
cdr3_dataframe_adaptive.fx <- function(datapath, chain, filelist, totalinframe) {
  
  filelist <- filelist[grepl(chain, filelist)]
  i <- 1
  for (f in filelist) {
    mixcrfle <- read.table(paste(datapath, f, sep = ""), 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                           na.strings = c("", "NA"))
    if (i == 1) {
      compldfle <- mixcrfle
      compldfle <- cbind(cloneno = row.names(compldfle), 
                         filename = f, compldfle)
      i <- i + 1
    }
    else {
      compldfle1 <- mixcrfle
      compldfle1 <- cbind(cloneno = row.names(compldfle1), 
                          filename = f, compldfle1)
      compldfle <- rbind(compldfle, compldfle1)
      rm(compldfle1)
    }
  }
  myfiles <- unique(as.character(compldfle$filename))
  message("my files:")
  print(myfiles)
  message("Total recovered clonotypes:")
  print(length(compldfle$aminoAcid))
  message("Total out-of-frame clonotypes:")
  print(length(compldfle$aminoAcid[compldfle$aminoAcid == ""]))
  compldfle$filename <- as.character(compldfle$filename)
  compldfle$samplename <- gsub(".*.CLONES_", "", compldfle$filename)
  compldfle$samplename <- gsub(chain, "", compldfle$samplename)
  compldfle_clean <- compldfle[compldfle$aminoAcid != "", ]
  compldfle_clean$cloneFraction <- NA
  for (f in myfiles) {
    compldfle_clean$cloneFraction[compldfle_clean$filename == 
                                    f] <- compldfle_clean$cloneCount[compldfle_clean$filename == 
                                                                       f]/sum(compldfle_clean$cloneCount[compldfle_clean$filename == 
                                                                                                           f])
  }
  message("Total productive clonotypes:")
  print(length(compldfle_clean$aminoAcid))
  if (totalinframe == "inframe") {
    message("Output contains in_frame clonotypes only")
    return(compldfle_clean)
  }
  if (totalinframe == "total") {
    message("Output contains all clonotypes")
    return(compldfle)
  }
}


# Track one specific nt seq using adaptive data
track_Aclone_adaptive.fx <- function(compldfle, plotpath, countfrac, clnefrc, clnseq){
  
  message("list of samples to track clones: ")
  mysamples <- unique(compldfle$samplename)
  print(mysamples)
  
  # Subset df
  CDR3_fraction <- compldfle[, c("samplename","nucleotide","cloneFraction", "cloneCount")]
  
  # Subset to include only clonotypes with more than specified clonal fraction
  CDR3_fraction <- CDR3_fraction[CDR3_fraction$cloneFraction > clnefrc,]
  
  #Assign colors to the specific clone
  myclone <- clnseq
  notrecurring <- CDR3_fraction$nucleotide[!CDR3_fraction$nucleotide %in% myclone]
  
  cloneColor <- distinctColorPalette(length(myclone))
  myColors <- c(cloneColor, rep("white",length(notrecurring)))
  names(myColors) <- c(myclone, notrecurring)
  
  # Generate a row for each sample that doesnot have jurkat clonotype
  ## This ensures alluvia are colored
  
  tmp <- CDR3_fraction[CDR3_fraction$nucleotide == myclone,]
  nonexisting <- mysamples[!mysamples %in% tmp$samplename]
  if(length(nonexisting) > 0){
    newentries <- data.frame("samplename" = nonexisting, "nucleotide" = myclone,
                             "cloneFraction" = 0, "cloneCount" = 0)
    CDR3_fraction <- rbind(CDR3_fraction, newentries)
  }
  
  p <-  ggplot(CDR3_fraction, aes(x = samplename,
                                  y = eval(as.name(countfrac)),
                                  fill = nucleotide,
                                  stratum = nucleotide,
                                  alluvium = nucleotide,
                                  label = nucleotide))
  
  myp <- p + geom_alluvium(decreasing = FALSE) +
    geom_stratum(decreasing = FALSE, stat = "alluvium") +
    scale_fill_manual(breaks = names(myColors[myColors != "white"]),
                      values = myColors) +
    theme(axis.title.y = element_text(size = 50),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.position = "none",
          plot.margin = unit(c(0.2,0,0,0),"cm")) +
    labs(y = countfrac)
  
  pdf(paste0(plotpath, "trackAclone_Adaptive_", mysamples[1], countfrac, ".pdf"),
      width = 30,
      height = 20,
      useDingbats = FALSE,
      onefile = FALSE)
  print(myp)
  dev.off()
  
}
