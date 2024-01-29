# round and format pvalues super large ones with sci notation
# 
round_and_format <- function(x, threshold = 0.001) {
  # Round to the second significant digit
  magnitude <- 10 ^ (floor(log10(abs(x))))
  rounded <- round(x / magnitude, 1) * magnitude
  
  # Format based on threshold
  if (rounded < threshold) {
    return(format(rounded, scientific = TRUE, digits = 2)) #two digits for scientific notation
  } else {
    return(format(rounded, scientific = FALSE, digits = 1)) # one digit for non-scientific notation
  }
}



#from SO: https://stackoverflow.com/questions/27627798/write-a-data-frame-to-a-xls-file-with-a-title

text_matrix <- function(dat, table_title) {
  rbind(c(table_title, rep('', ncol(dat)-1)), # title
        rep('', ncol(dat)), # blank spacer row
        names(dat), # column names
        unname(sapply(dat, as.character))) # data
}


# Add jurkat to make test dataset
addjurkat.fx <- function(f1, Jurkat, percincrement, outpath){
  jurkatclone <- Jurkat[1,]
  print(summary(f1$cloneCount))
  print(sum(f1$cloneCount))
  increment <- sum(f1$cloneCount) * percincrement
  myseq <- seq(0, sum(f1$cloneCount), increment)
  
  #get random rows
  set.seed(777)
  samprows <- sample(1:nrow(f1))     
  for(i in myseq[myseq != 0]){
    message(i)
    jurkatclone$cloneCount <- i
    myperc <- i / sum(f1$cloneCount)
    # random rows that cumsum is close to i
    removerows <- samprows[which(cumsum(f1[samprows, "cloneCount"]) <= i)] 
    print(sum(f1$cloneCount[removerows]))
    # bind jurkat clone to f1
    mycloneset <- rbind(f1[-removerows,], jurkatclone)
    #recalculate clonefraction
    mycloneset$cloneFraction <- mycloneset$cloneCount/sum(mycloneset$cloneCount)
    write.table(mycloneset, file = paste0(outpath, "CLONES_TRBaddjurkat_", i,"_",myperc, ".txt"), quote = F, sep = "\t")   
  }
}


#plot histograms
histp <- function (df, var, bin) {
  myp <- ggplot(data = df, aes(x = eval(parse(text = var)))) + 
    geom_histogram(bins = bin) + myplot + myaxis + labs(x = var)
  return(myp)
}

#plot histograms and density
histdenp <- function (df, var, bin) {
  myp <- ggplot(data = df, aes(x = eval(parse(text = var)))) + 
    geom_histogram(aes(y = after_stat(density)), fill = "white", color = "black" , binwidth = bin) +
    geom_density(color = "red", size = 1) + myplot + myaxis + labs(x = var)
  return(myp)
}


#to align plots (from stackoverflow)
# Function to align plots (from stackoverflow) 
align_plots1 <- function (...) {
  pl <- list(...)
  stopifnot(do.call(all, lapply(pl, inherits, "gg")))
  gl <- lapply(pl, ggplotGrob)
  bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
  combined <- Reduce(bind2, gl[-1], gl[[1]])
  wl <- lapply(gl, "[[", "widths")
  combined$widths <- do.call(grid::unit.pmax, wl)
  grid::grid.newpage()
  grid::grid.draw(combined)
}




#mclapply to print warnings and errors from SO: https://stackoverflow.com/questions/21486658/warnings-suppressed-with-mclapply-in-r
safe_mclapply <- function(X, FUN, mc.cores, stop.on.error=T, ...){
  fun <- function(x){
    res_inner <- tryCatch({
      withCallingHandlers(
        expr = {
          FUN(x, ...)
        }, 
        warning = function(e) {
          message_parallel(trimws(paste0("WARNING [element ", x,"]: ", e)))
          # this line is required to continue FUN execution after the warning
          invokeRestart("muffleWarning")
        },
        error = function(e) {
          message_parallel(trimws(paste0("ERROR [element ", x,"]: ", e)))
        }
      )},
      error = function(e){
        # error is returned gracefully; other results of this core won't be affected
        return(e)
      }
    )
    return(res_inner)
  }
  
  res <- mclapply(X, fun, mc.cores=mc.cores, mc.preschedule = FALSE)
  failed <- sapply(res, inherits, what = "error")
  if (any(failed == T)){
    error_indices <- paste0(which(failed == T), collapse=", ")
    error_traces <- paste0(lapply(res[which(failed == T)], function(x) x$message), collapse="\n\n")
    error_message <- sprintf("Elements with following indices failed with an error: %s. Error messages: \n\n%s", 
                             error_indices,
                             error_traces)
    if (stop.on.error)
      stop(error_message)
    else
      warning(error_message, "\n\n### Errors will be ignored ###")
  }
  return(res[!failed])
}

#' Function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}


# Add a specific clone (eg jurkat) incrementally to generate some test dataset
addjurkat.fx <- function(f1, Jurkat, percincrement, outpath){
  jurkatclone <- Jurkat[1,]
  print(summary(f1$cloneCount))
  print(sum(f1$cloneCount))
  increment <- sum(f1$cloneCount) * percincrement
  myseq <- seq(0, sum(f1$cloneCount), increment)
  
  #get random rows
  set.seed(777)
  samprows <- sample(1:nrow(f1))     
  for(i in myseq[myseq != 0]){
    message(i)
    jurkatclone$cloneCount <- i
    myperc <- i / sum(f1$cloneCount)
    # random rows that cumsum is close to i
    removerows <- samprows[which(cumsum(f1[samprows, "cloneCount"]) <= i)] 
    print(sum(f1$cloneCount[removerows]))
    # bind jurkat clone to f1
    mycloneset <- rbind(f1[-removerows,], jurkatclone)
    #recalculate clonefraction
    mycloneset$cloneFraction <- mycloneset$cloneCount/sum(mycloneset$cloneCount)
    #order clonefraction
    mycloneset <- mycloneset[ order(mycloneset$cloneCount, decreasing = T),]        
    # make clone count integer
    mycloneset$cloneCount <- as.integer(mycloneset$cloneCount) 
    write.table(mycloneset, file = paste0(outpath, "CLONES_TRBaddjurkat_", i,"_",myperc, ".txt"), quote = F, sep = "\t")   
  }
}

sampletags_columns <- function(orig_df, grepvars) {
  orig_df$index <- 1:nrow(orig_df) # add index column
  orig_df$sample_tags <- paste0(orig_df$index,",",orig_df$sample_tags) # add index to sample_tags
  splitsampletags <- strsplit(orig_df$sample_tags, split = ",") # split sample_tags by comma into a list
  # for each sample_tag, extract the variables in grepvars as list
  mydf <- lapply(splitsampletags, function(sampletag){ 
    y <- unlist(sampletag)
    y <- trimws(y)
    indx <- y[1] # first element is index
    # for each variable in grepvars, extract the value if it exists, if not add NA
    myvars <- lapply(grepvars, function(myvar){
      ifelse(sum(grepl(myvar, y)) == 1, y[grepl(myvar, y)], NA) })
    myvarsdf <- as.data.frame(myvars)
    vardf <- cbind.data.frame(indx, myvarsdf)
    colnames(vardf) <- c("index", grepvars) # rename columns
    return(vardf)
  })
  return(do.call(rbind, mydf)) # return a data frame
}



toString_onefle.fx <- function(df, tostringvar){
  setDT(df)
  ab <- df[, .(tostringvar = toString( eval(parse( text = tostringvar))), #bind subjects
               count = sum(count)),  #get sum of counts
           by = c("CDR3b", "TRBV", "TRBJ")] #get duplicates sequences with the same cdr3 + TRBV + TRBJ
  df_ab <- merge(df, ab, by = c("CDR3b","TRBV","TRBJ")) # merge together
  df_ab_dedup <- dplyr::distinct(df_ab, CDR3b, TRBV, TRBJ, tostringvar, .keep_all= TRUE)
  
  return(df_ab_dedup[, c("CDR3b", "TRBV", "TRBJ", "samplename", 
                         "count.y", "clonefraction", "tostringvar", "file")])
}

# spiderplot function
calculate_delta.fx <- function(df, var1, var2) {
  # make a tbale patient x var1
  mytab <- table(df$Patient, df[[var1]])
  # min two samples per patient
  mytab <- mytab[rowSums(mytab == 1) > 1, ]
  # select those samples with a baseline
  baseline_patients <- rownames(mytab)[mytab[, 1] == 1]
  df1 <- df[df$Patient %in% baseline_patients, ]
  result <- df1 %>%
    group_by(Patient) %>%
    mutate(Difference = eval(parse(text = var2)) - eval(parse(text = var2))[eval(parse(text = var1)) == "X01"])
  return(result)
}


delta_basespiderplot.fx <- function(df_diff, var1, clrby, colpal) {
  # df_diff from calculate_delta.fx
  p0 <- ggplot(
    df_diff,
    aes(x = eval(parse(text = var1)), y = Difference)
  ) +
    geom_point(aes(color = eval(parse(text = clrby))), cex = 2) +
    geom_line(aes(group = Patient, color = eval(parse(text = clrby)))) +
    scale_color_manual(values = colpal) +
    myplot +
    myaxis +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  return(p0)
}




