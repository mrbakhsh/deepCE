


  #' getPPI
  #' @title Generating Protein-Protein Interaction (PPI) Vector
  #' @param data A Co-elution data matrix with proteins in rows and
  #' fractions in columns.
  #' @param similarity_calculation If TRUE, this function first computes the
  #' similarity of protein pairs based on co-elution pattern via one of the
  #' user-defined metrics. It then removes protein pairs below 0.5
  #' similarity score
  #' @param metric The measure of association including:
  #' c(pearson", "spearman", "kendall", "bicor", "cosine", "jaccard",
  #' "euclidean", "manhattan", "weighted_rank", "mutual_information",
  #' "zero_count_bayesian", "weighted_cross_correlation", or
  #' "pearson_plusNoise"); Defaults to pearson.
  #' See \code{\link{similarity_score}}.
  #' @return A concatenated matrix.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function concatenates co-elution profiles for each
  #' possible protein pairs to form a master matrix for input into the learner.
  #' @importFrom dplyr bind_cols
  #' @export
  #' @examples
  #' data("HelaCE")
  #' # select subset of a data
  #' m <- HelaCE[1:10, 1:10]
  #' conc_dat <- getPPI(m, similarity_calculation = FALSE)



  getPPI <-
    function(data,
             similarity_calculation = FALSE,
             metric = "pearson"){


      if (!is.matrix(data)) {
        data <- as.matrix(data)
      }

      if(is.character(data) == TRUE){
        stop("matrix must include numerical variables")
      }

      if(is.null(row.names(data))){
        stop("Please specify the row.names")
      }

      rowname <- NULL
      value <- NULL

      data <- filter_ConsecutivePep(data, 2)


      if(similarity_calculation){

        data_s <- data[sort(rownames(data)), ]
        m <-  similarity_score(data_s, metric)
        m[lower.tri(m, diag=TRUE)] <- NA
        s <-
          m %>% as.data.frame %>% rownames_to_column() %>%
          pivot_longer(-rowname) %>%
          na.omit() %>%
          filter(value > 0.5) %>%
          unite(PPI, c("rowname", "name"), sep = "~", remove = FALSE)
        colnames(s)[2:3] <- c("p1", "p2")
        p1 <-
          s %>% select("p1")

        p2 <-
          s %>% select("p2")

        dat <- as.data.frame(data_s)
        dat <- rownames_to_column(dat, "p1")

        out_p1<-plyr::join(p1, dat, by = "p1")
        out_p1 <- out_p1[, -1]
        colnames(dat)[1] <- "p2"
        out_p2<-plyr::join(p2, dat, by = "p2")
        out_p2 <- out_p2[, -1]

        co_elutdf <- bind_cols(out_p1,out_p2)
        co_elutdf <- as.matrix(co_elutdf)
        row.names(co_elutdf) <- s$PPI
        message("Number of retained PPIs after correlation removal:",
                nrow(co_elutdf))

      }
      if(similarity_calculation ==  FALSE) {
      r <- data[sort(rownames(data)), ]
      PPI <-
        as.data.frame(comboGeneral(rownames(r), 2))
      colnames(PPI) <- c("p1", "p2")
      pair <- unite(PPI, PPI, c(p1,p2), sep = "~")

      p1 <-
        PPI %>% select("p1")

      p2 <-
        PPI %>% select("p2")


      data <- as.data.frame(data)
      data <- rownames_to_column(data, "p1")

      out_p1<-plyr::join(p1, data)
      out_p1 <- out_p1[, -1]
      colnames(data)[1] <- "p2"
      out_p2<-plyr::join(p2, data)
      out_p2 <- out_p2[, -1]


      co_elutdf <- bind_cols(out_p1,out_p2)
      co_elutdf <- as.matrix(co_elutdf)
      row.names(co_elutdf) <- pair$PPI
      }

      return(co_elutdf)
    }
