# 
# S plot similar to Lawrence paper
#
median.cohorts.fx <- function(df, myvar, mygroups){
  df.mygroups <- cbind(unique(df[[mygroups]]), NA)
  colnames(df.mygroups) <- c("group", "median") 
  df.mygroups <- as.data.frame(df.mygroups) 
  df.mygroups$median <- as.numeric(df.mygroups$median)
  for(i in 1:nrow(df.mygroups)){
    
    df.mygroups$median[i]<-median(df[[myvar]][df[[mygroups]] == df.mygroups$group[i]])
  }
  df.mygroups <- df.mygroups[order(df.mygroups$median, decreasing = F),]
  return(df.mygroups)
}


sort.df.fx <- function(df, median_df, myvar, mygroups){
  disease.width <- (nrow(df)/nrow(median_df)) 
  sorted.df <- df[0,]
  start = 0
  for(i in 1:(nrow(median_df))){
    tmp <- df[df[[mygroups]]==median_df$group[i],]
    tmp <- tmp[order(tmp[[myvar]]),]
    tmp <- tmp[!is.na(tmp[[myvar]]),] 
    #create range of x values to squeeze dots into equal widths of the plot for each Disease regardless of the number of samples
    div <- disease.width/nrow(tmp)
    #If there is only one sample, put the dot in the middle of the alloted space
    if(dim(tmp)[1]==1){
      tmp$Xpos<-start+(disease.width/2)
    } 
    else tmp$Xpos <- seq(from = start, 
                         to = start+disease.width, 
                         by = div)[-1]
    sorted.df <- rbind(sorted.df, tmp)  
    median_df$Median.start[i] <- tmp$Xpos[1]
    median_df$Median.stop[i] <- tmp$Xpos[nrow(tmp)]
    median_df$N[i]<-nrow(tmp)
    start <- start+disease.width+30
    
  }
  median_df$medianloc <- median_df$Median.start+
    ((median_df$Median.stop-median_df$Median.start)/2)
  sorted.df$group <- factor(sorted.df[[mygroups]],
                            levels = median_df$group)
  return(list(sorted.df,median_df))
}


Splot.fx <- function(list.sorted_df.median, myvar, colby, colpal, plottitle){       
  sorted_df <- as.data.frame(list.sorted_df.median[[1]])
  median_df <- as.data.frame(list.sorted_df.median[[2]])
  disease.width <- (nrow(sorted_df)/nrow(median_df)) 
  
  
  Splot <- ggplot() +
    geom_point(data = sorted_df, 
               aes(x = Xpos ,y = eval(parse(text = myvar)),
                   color = eval(parse(text = colby))), 
               size = 5) +
    geom_crossbar(data = median_df, 
                  aes(x =medianloc, y = median,
                      ymin = median, ymax = median),
                  width = disease.width) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 22, color = "black"),
          axis.title = element_text(size = 22), 
          plot.title = element_text(size=22, hjust = 0.5),
          legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",
                                          colour = NA),
          panel.border=element_blank(),
          plot.margin = unit(c(1.2,1,0,1),"cm")) +
    scale_color_manual(values = colpal, guide = FALSE) +
    scale_x_continuous(breaks = seq((disease.width)/2,max(sorted_df$Xpos),
                                    disease.width+30),
                       labels = median_df$group,
                       expand = c(0,20)) + 
    labs(y = myvar, title = plottitle) 
  return(Splot)  
}

# plot PCA related to Fig3
#
dimplot_pt <- function(df, pt, myshape, myalpha) {
  p1 <- ggplot(data = df, aes(x = Dim1, y = Dim2, shape = eval(parse(text = myshape)) )) +
    geom_point(aes(color = cancergroup, alpha = eval(parse(text = myalpha))), size = 3) +
    geom_path(
      data = df[df$Patient == pt, ], aes(group = Patient), # color = c("red", "black"),
      arrow = arrow(length = unit(0.30, "cm"), ends = "last", type = "closed")
    ) +
    myplot +
    myaxis +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    scale_color_manual(values = group_col) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggtitle(pt)
  return(p1)
}





