    .findMissingValues <- function(data,
                                  bound_left = 2,
                                  bound_right = 2){
      intMat <- data
      consider_borders = TRUE
      #Pad matrix for border neighbors
      padVal <- ifelse(consider_borders, NA, 0)
      padLeft <- matrix(padVal, nrow= nrow(intMat), ncol = bound_right)
      padRight <- matrix(padVal, nrow= nrow(intMat), ncol = bound_left)
      intMatPad <- cbind(padLeft, intMat, padLeft)
      # Find missing values
      zeroes <- which(intMatPad == 0, arr.ind = TRUE)
      zeroes <- zeroes[(zeroes[,2] > bound_left) &
                         (zeroes[,2] <= (ncol(intMatPad) - bound_right)),]
      zeroMat <- matrix(0, nrow=nrow(intMat),
                        ncol = ncol(intMat)+bound_left+bound_right)
      neighbors <- c(-bound_left:-1, 1:bound_right)
      z <- zeroes
      for(i in neighbors){
        z[,2] <- zeroes[,2] + i
        zeroMat[z] <- zeroMat[z] +1
      }

      intMatPad[,c(1:bound_left,
                   (ncol(intMatPad)-bound_right+1):ncol(intMatPad))] <- 999
      missingVals <- which(intMatPad == 0 & zeroMat == 0, arr.ind = TRUE)
      missingVals[,2] <- missingVals[,2] - bound_left
      data[missingVals] <- NA
      return(data)
    }

    .imputeMissingVals <- function(data){

      intMat <- data
      naIndx <- which(is.na(intMat), arr.ind = TRUE)

      intMatImp <- apply(intMat, 1, function(tr){
        naIdx <- is.na(tr)
        n <- length(tr)

        if (!any(naIdx)) {
          return(tr)
        }

        allindx <- 1:n
        indx <- allindx[!naIdx]

        interp <- spline(indx, tr[indx], n = n)$y

        # Calculate border fractions (linear interpolation of 2
        # adjacent fractions)
        if(naIdx[1]){
          interp[1] <- 2 * interp[2] -interp[3]
        }
        if(naIdx[n]){
          interp[n] <- 2 * interp[n-1] -interp[n-2]
        }

        return(interp)
      })

      intMatImp <- t(intMatImp)
      intMatImp[intMatImp < 0] <- 0
      data[naIndx] <- round(intMatImp[naIndx],digits=2)
      imputed_data <- data
      return(imputed_data)
    }


    .Plot_NAs <- function(data_na){

      Fraction <- NA
      NAs <- NA

      intMat <- data_na
      ## NA distribution
      #per fraction
      nas <- apply(intMat, 2, function(x) length(which(is.na(x))))

      p <-
        ggplot(data.table(Fraction = 1:length(nas), NAs = nas)) +
        geom_bar(aes(x = Fraction, y = NAs),stat = "identity") +
        ylab("Missing values") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        theme(text = element_text(size = 12))   +
        ggtitle("Missing Values per fraction") +
        theme_classic() +
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

    .Plot_imputations <- function(data_na, data_imputed){

      . <- NA
      meanV <- NA
      type <- NA
      value <- NA
      key <- NA

    prior_imp <-
      as.data.frame(data_na) %>%
      gather(key, value, 1:ncol(.)) %>%
      group_by(key) %>% summarise(meanV = mean(value, na.rm = TRUE)) %>%
      mutate(type = "prior imputation")

    after_imp <-
      as.data.frame(data_imputed) %>%
      gather(key, value, 1:ncol(.) )%>%
      group_by(key) %>% summarise(meanV = mean(value, na.rm = TRUE)) %>%
      mutate(type = "after imputation")

     df_plot <- rbind(prior_imp, after_imp)
     p <-
       ggplot(df_plot, aes( x = meanV, color=type)) +
      geom_density(size = 1) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      ylab("Density") +
      xlab("Average intesity value per fraction") +
      ggtitle("Missing Value Imputation") +
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


  #' impute_MissingData
  #' @title Detect and Impute Missing Values
  #' @param data A data matrix with rows
  #' including proteins and fractions along the columns.
  #' @param bound_left Numeric integer, the minimum number of non-zero
  #' values to the left of a missing value to be replaced with \code{NA};
  #' defaults to 2.
  #' @param bound_right Numeric integer, the minimum number of non-zero
  #' values to the right of a missing value to be replaced with \code{NA};
  #' defaults to 2.
  #' @param plotImputationSummary Logical value, indicating whether to  plot
  #' the imputation summary; defaults to FALSE.
  #' @param filename Character string, indicating the output file name as an
  #' pdf object.
  #' @return Imputed matrix.
  #' @references
  #' Bludau,I. et al. (2020) Complex-centric proteome profiling by SEC-SWATH-MS
  #' for the parallel detection of hundreds of protein complexes.
  #' Nat. Protoc., 15, 2341-2386.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom data.table data.table
  #' @importFrom ggplot2 geom_bar
  #' @importFrom ggplot2 ggtitle
  #' @importFrom ggplot2 theme_classic
  #' @importFrom ggplot2 geom_density
  #' @importFrom ggplot2 scale_x_continuous
  #' @importFrom tidyr gather
  #' @importFrom stats spline
  #' @description The algorithm first identifies 0 values on the data set as
  #' missing values if they fulfill the rule:
  #' \code{(1)*bound_left - 0 - (1)*bound_right} i.e. a zero value has to have
  #' at least a specified number of non-zero neighbors to be classified
  #' as \code{NA}. NA values are then linearly extrapolated from the
  #' neighboring  values via interpolation. Any imputed value below 0
  #' is set to 0.
  #' @export
  #' @examples
  #' # Load the data
  #' data("HelaCE")
  #' # Impute the missing values
  #' Imputed_data <-
  #' impute_MissingData(HelaCE,bound_left = 2,bound_right = 2)



  impute_MissingData <-
    function(data,
             bound_left = 2,
             bound_right = 2,
             plotImputationSummary = FALSE,
             filename = "plots.pdf") {

    if (!is.matrix(data)) {
      data <- as.matrix(data)
    }


    if(is.character(data) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(data))){
      stop("Please specify the row.names")
    }




    data_na <- .findMissingValues(data)

    data_imputed <- .imputeMissingVals(data_na)

    if (plotImputationSummary) {

      pdf(filename)

      na_plot <-
        .Plot_NAs(data_na)
      imputation_summary <-
        .Plot_imputations(data_na, data_imputed)

    # Make plots.
    plot_list <- list()
    plot_list$na_plot <- na_plot
    plot_list$summary <- imputation_summary

    for (i in seq_along(plot_list)) {
        print(plot_list[[i]])
    }
      dev.off()

    }

    return(data_imputed)
  }
