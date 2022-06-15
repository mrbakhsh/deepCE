  .Plot_RemovingSinglePeaks <- function(data, dataConsPeaks_m){

    . <- NA
    meanV <- NA
    type <- NA
    value <- NA
    key <- NA

    prior_deletion <-
      as.data.frame(data) %>%
      gather(key, value, 1:ncol(.)) %>%
      group_by(key) %>% summarise(meanV = mean(value, na.rm = TRUE)) %>%
      mutate(type = "original")

    after_deletion <-
      as.data.frame(dataConsPeaks_m) %>%
      gather(key, value, 1:ncol(.) )%>%
      group_by(key) %>% summarise(meanV = mean(value, na.rm = TRUE)) %>%
      mutate(type = "removed single peaks")

    df_plot <- rbind(prior_deletion, after_deletion)
    p <-
      ggplot(df_plot, aes( x = meanV, color=type)) +
      geom_density(size = 1) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      ylab("Density") +
      xlab("Average intesity value per fraction") +
      ggtitle("Remove single peaks") +
      theme_classic() +
      theme(text = element_text(size = 12))   +
      theme( axis.line = element_line(size = 0.75))  +
      theme(axis.text = element_text(
        color = "black",
        size = 12,colour = "black")) +
      theme(axis.ticks = element_line(
        colour = "black",
        size = 0.75
      )) +
      theme(axis.ticks.length=unit(.3, "cm"))
  }


  #' RemoveSinglePeak
  #' @title Remove Single Peak
  #' @param data A data matrix with rows including proteins and fractions
  #' along the columns.
  #' @param plot_RemovalSinglePeaks Logical value, indicating whether to  plot
  #' the distribution after peak removal; Defaults to FALSE.
  #' @param filename Character string, indicating the output file name as an
  #' pdf object.
  #' @return Filtered matrix.
  #' @references
  #' Schlossarek,D. et al. (2021) PROMISed: A novel web-based tool to
  #' facilitate analysis and visualization of the molecular interaction
  #' networks from co-fractionation mass spectrometry (CF-MS) experiments.
  #' Comput. Struct. Biotechnol.
  #' @description This function replaces single values surrounded by zeroes
  #' with zero, eliminating data-noise “peaks” that only span one fraction.
  #' @export
  #' @examples
  #' # load the co-eluton table
  #' data("HelaCE")
  #' # select subset of a data
  #' m <- HelaCE[1:10,1:10]
  #' filt_data <- RemoveSinglePeak(m)

  RemoveSinglePeak <- function(data,
                               plot_RemovalSinglePeaks = FALSE,
                               filename =  "plots.pdf"){

     if (is.matrix(data)) {
      data <- as.data.frame(data)
    }

    if(is.character(data) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(data))){
      stop("Please specify the row.names")
    }

   z <- data

    dataConsPeaks <- as.data.frame(NULL)
    w <- 1
    while (w <= nrow(z)){
      tmpNos <- z[w,] # copy intensities from row w of data x
      tmpCount <- which(tmpNos>0) # check which "cells" are higher than 0
      tmp<-c() # create empty vector
      tmpr<-c() # create empty vector

      for (p in 1:length(tmpCount)){
        tmp[p]<- tmpCount[p+1]-(tmpCount[p])
      }

      for (q in 1:length(tmpCount)){
        tmpr[q]<- tmpCount[q]-(tmpCount[q-1])
      }

      tmpb <- NULL
      for(i in 1:length(tmp)){
        tmpb[i] <- tmp[i] == 1 | tmpr[i] == 1
      }

      tmpb[is.na(tmpb)] <- FALSE

      tmpCountReal <- tmpCount[tmpb]

      tmpm <- as.data.frame(matrix(0, nrow = nrow(tmpNos), ncol = ncol(tmpNos)))
      x <- 1
      y <- 1

      while(x <= length(tmpNos)){
        while(y <= length(tmpCountReal)){
          if(x == tmpCountReal[y]){
            tmpm[,x] <- tmpNos[x]
            y <- y+1
          }else{
            y <- y+1
          }
        }
        y <- 1
        x <- x+1
      }

      colnames(tmpm) <- colnames(tmpNos)
      rownames(tmpm) <- rownames(tmpNos)

      dataConsPeaks <- rbind(dataConsPeaks, tmpm)
      w <- w+1
    }

    dataConsPeaks_m <-
      as.matrix(dataConsPeaks)

    if (plot_RemovalSinglePeaks) {

      pdf(filename)

      remoPeaks_plot <-
        .Plot_RemovingSinglePeaks(data,dataConsPeaks_m)

      print(remoPeaks_plot)
      dev.off()
      }

    return(dataConsPeaks_m)
  }
