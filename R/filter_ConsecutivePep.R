  #' filter_ConsecutivePep
  #' @title Filter Peptides by Consecutive Peptide Detection
  #' @param data A data matrix with rows including proteins and fractions
  #' along the columns.
  #' @param min_stretch_length Numeric integer, the minimal length a stretch of
  #' continuous identifications has to have in order not to be removed.
  #' Defaults to 2.
  #' @return Filtered matrix
  #' @references
  #' Bludau,I. et al. (2020) Complex-centric proteome profiling by SEC-SWATH-MS
  #' for the parallel detection of hundreds of protein complexes.
  #' Nat. Protoc., 15, 2341-2386.
  #' @description This function removes peptides/proteins that have never been
  #' detected in more than N consecutive fractions (here N = 2).
  #' @export
  #' @examples
  #' # Load the data
  #' data("HelaCE")
  #' # Remover proteins that detected in more than 10 consecutive fractions
  #' filt_dat <- filter_ConsecutivePep(HelaCE,10)

  filter_ConsecutivePep <-
    function(data, min_stretch_length=2) {

      if (!is.matrix(data)) {
        data <- as.matrix(data)
      }

      if(is.character(data) == TRUE){
        stop("matrix must include numerical variables")
      }

      if(is.null(row.names(data))){
        stop("Please specify the row.names")
      }

      dummy <- NULL


      ## Get traces from container
      id <- rownames(data)
      intensity <- data
      ## Add 0-column
      intensity <- cbind(intensity, dummy=rep(0, nrow(intensity)))

      ## Define n as the number of columns
      n <- ncol(intensity)
      ## Set count variable to one
      tmp <- 1
      ## Going through all rows, for all do:
      for (x in seq_len(nrow(intensity))) {
        tmp <- 1
        for (i in seq_len(n)) {
          if (intensity[x,i] == 0) {
            tmp <- 1
          } else {
            if (intensity[x, i+1] == 0) {
              if (tmp <= (min_stretch_length-1)){
                for (j in 0:tmp-1) {
                  intensity[x, i-j] <- 0
                }
                tmp <- 1
              }
            } else {
              tmp <- tmp + 1
            }
          }
        }
      }

      # Remove dummy column
      intensity <- subset(intensity, select =-dummy)
      data_filtered <-
        intensity[rowSums(intensity) != 0,]
      return(data_filtered)

    }











