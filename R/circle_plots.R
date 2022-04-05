circlepack.reads.fx <- function(inputfile, sample_id_DNA){
  
  sample_df <- inputfile[inputfile$samplename == sample_id_DNA,]
# type is gliph groups that I added to mixcr output file. This table gives
  type_tab <- as.data.frame(table(sample_df$Type), stringsAsFactors = F)

  if(nrow(type_tab) != 0){    
    cluster_tab <- type_tab[type_tab$Freq > 1,]
  } else {
    cluster_tab <- type_tab 
  }    
  
  # if there is a cluster, group cdr3s and use sum of reads
  if(nrow(cluster_tab) != 0){    
    cluster_tab$typereads <- NA
    for(j in 1:nrow(cluster_tab)){
      mytype <- cluster_tab$Var1[j]
      mytypereads <- sum(sample_df$cloneCount[which(sample_df$Type == mytype)])
      cluster_tab$typereads[j] <- mytypereads
    }
    #clean up cluster_tab
    cluster_tab <- cluster_tab[, c("Var1", "typereads")]
    colnames(cluster_tab) <- c("name", "size")
  }
  
  # Make edge df
  type_cdr3 <- sample_df[, c("Type", "nSeqCDR3")]
  colnames(type_cdr3) <- c("from", "to")
  #if no clusters, edge df is octamer_cdr3 df and all octamers will be replaced by sample_id
  if(nrow(cluster_tab) == 0) {
    myedge <- type_cdr3
    myedge[,"from"] <- unique(sample_df$samplename)
  }
  
  # if there is a cluster, make a list of sample_id and clusters
  if(nrow(cluster_tab) != 0){
    myclusters <- cbind.data.frame(NA,cluster_tab$name, stringsAsFactors = F)
    colnames(myclusters) <- c("from", "to")
    myclusters$from <- unique(sample_df$samplename)
    # replace types < 2 sequences with sample_id
    type_cdr3$from[!type_cdr3$from %in% cluster_tab$name] <- unique(sample_df$samplename)
    myedge <- rbind(myclusters, type_cdr3)   
  }
  
  # get cdr3 and reads    
  cdr3_freq <- sample_df[, c("nSeqCDR3", "cloneCount")]
  colnames(cdr3_freq) <- c("name", "size")
  
  #bind all and cleanup
  # get sample frequency        
  sample_tab <- as.data.frame(table(sample_df$samplename), stringsAsFactors = F)
  colnames(sample_tab) <- c("name", "size")
  
  #bind all and cleanup
  myvertex <- rbind(sample_tab, cdr3_freq)
  
  # if there is a cluster, include it to vertex
  if(nrow(cluster_tab) != 0){
    myvertex <- rbind(myvertex, cluster_tab)
  }
  
  # first row is NA so remove it
  myvertex <- myvertex[!is.na(myvertex$name),]
  myvertex$size <- as.numeric(myvertex$size)
  
  # Make a type variable for colors
  myvertex$type <- NA
  myvertex$type[1] <- "samplename"
  myvertex$type[myvertex$name %in% cdr3_freq$name] <- "CDR3"
  myvertex$type[is.na(myvertex$type)] <- myvertex$name[is.na(myvertex$type)] 
  
  
  myvertex$group <- NA
  myvertex$group[myvertex$type == "samplename"] <- "samplename"
  myvertex$group[myvertex$type == "CDR3"] <- "CDR3"
  myvertex$group[!myvertex$type %in% c("CDR3", "samplename")] <- "Type"
  
  if(nrow(cluster_tab) != 0){
    myColors <- distinctColorPalette(nrow(cluster_tab))
    myColors <- c("white", # color samplename white
                  rep("#ededed",nrow(myvertex[myvertex$type == "CDR3",])), #color all cdr3 light gray
                  myColors)
    names(myColors) <- myvertex$type
  } else{
    myColors <- c("white", 
                  rep("#ededed",nrow(myvertex[myvertex$type == "CDR3",])))    
    names(myColors) <- myvertex$type
  }
  
  mygraph <- graph_from_data_frame(myedge, vertices = myvertex)  
  
  alphapal <- c("CDR3" = 1,"Type" = 0.4, "samplename" = 0)
  circlep <- ggraph(mygraph, layout = 'circlepack', weight = size) + 
    geom_node_circle(aes(fill = type, alpha = group)) +
    theme_void() + scale_fill_manual(values = myColors) +
    scale_alpha_manual(values = alphapal)
  
  
  pdf(file = paste0(plotpath,sample_id_DNA,"_gliph_circles.pdf"),
      width = 15, 
      height = 15,
      useDingbats = FALSE)
  print(circlep)
  dev.off()
}