  .conc <- function(prot1, prot2)
  {
    if (!is.matrix(prot1)) prot1 <- as.matrix(prot1)
    if (!is.matrix(prot2)) prot2 <- as.matrix(prot2)

    prot1Row <- nrow(prot1)
    prot2Row <- nrow(prot2)

    if (prot1Row != prot2Row) stop("Matrix row count must match")


    result <- cbind(prot1, prot2)


    return(result)
  }



  #' build_trainingData
  #' @title Construct Training Data
  #' @param data A co-elution data matrix with proteins in rows and
  #' fractions in columns.
  #' @param refcpx A list of protein complexes to create class labels for
  #' training set.
  #' @return A list containing training data, along with class labels.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom RcppAlgos comboGeneral
  #' @importFrom dplyr sample_n
  #' @description This function creates a training set with class labels,
  #' where "1" corresponds to positive interactions (i.e., interactions
  #' belongs to  the same complex), and "0" corresponds to negative
  #' interactions (i.e., random interactions). To avoid prediction biases,
  #' this function creates a class labels with the same number of pairs. The
  #' training set then can be used as input for model tuning and model training.
  #' @export
  #' @examples
  #' # Load the list of ground-truth complexes
  #' data("refcpx")
  #' # Load the co-elution data
  #' data("HelaCE")
  #' # Buid the training set
  #' t_data <- build_trainingData(HelaCE, refcpx)



  build_trainingData <-
    function(data, refcpx){


    if (!is.matrix(data)) {
      data <- as.matrix(data)
    }

    if(is.character(data) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(data))){
      stop("Please specify the row.names")
    }

    if(!is.list(refcpx)){
      stop("Refcpx Must Be List")
    }

    # create a reference set
    r <- data[sort(rownames(data)), ]
    PPI <-
      as.data.frame(comboGeneral(rownames(r), 2))
    colnames(PPI) <- c("p1", "p2")

    class_labels <-
      generate_refInt(PPI,refcpx)

    # 1:1 positive/negative ratio
    pos <- filter(class_labels, labels == 1)
    neg <-
      class_labels %>%
      filter(labels == 0) %>% sample_n(nrow(pos))
    class_labels <- rbind(pos,neg)



    # create a training set #
    rnames <- row.names(data)
    p1 <- matrix(NA, nrow = nrow(class_labels), ncol = ncol(data))
    for (i in seq_len(nrow(class_labels)))
      p1[i, ] <- data[which(class_labels$p1[i] == rnames), ]

    # second profile
    p2 <- matrix(NA, nrow = nrow(class_labels), ncol = ncol(data))
    for (i in seq_len(nrow(class_labels)))
      p2[i, ] <- data[which(class_labels$p2[i] == rnames), ]


    # create training set
    train_d <- .conc(p1,p2)
    # create a numeric label set
    train_l <- class_labels$labels


    output_list <- list()
    output_list[["train_d"]] <- train_d
    output_list[["train_l"]] <- train_l

    return(output_list)


  }
