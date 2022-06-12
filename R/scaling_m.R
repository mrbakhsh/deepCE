  #' scaling_m
  #' @title Row-wise Normalization
  #' @param data A data matrix with rows
  #' including proteins and fractions along the columns.
  #' @return Scaled data matrix.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function performs row-wise normalization by dividing
  #' the ion intensity or spectral count of the particular protein by the
  #' max of ion intensities or spectra counts across fractions.
  #' @export
  #' @examples
  #' data("HelaCE")
  #' n_data <- scaling_m(HelaCE)



  scaling_m <- function(data) {

    if (!is.matrix(data)) {
      data <- as.matrix(data)
    }

    if(is.character(data) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(data))){
      stop("Please specify the row.names")

    }
    scaleM <- apply(data, 1,
                      function(x) (x)/max(x))
    scaleM <- t(scaleM)


    return(scaleM)

  }
