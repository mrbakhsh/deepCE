
  .cm_f <- function(preds, real){
    cm <-
      table(preds, real)

    TP <- as.double(cm[1, 1])
    TN <- as.double(cm[2, 2])
    FP <- as.double(cm[1, 2])
    FN <- as.double(cm[2, 1])

    ACC <- (TP + TN) / (TP + TN + FP + FN)
    SE <- TP / (TP + FN)
    SP <- TN / (FP + TN)
    PPV <- TP / (TP + FP)
    F1 <- 2 * TP / (2 * TP + TP + FN)
    MCC <- (TP * TN - FP * FN) /
      sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    cm_result <- list()

    cm_result[["cm"]] <- cm
    cm_result[["ACC"]] <- ACC
    cm_result[["SE"]] <- SE
    cm_result[["SP"]] <- SP
    cm_result[["PPV"]] <- PPV
    cm_result[["F1"]] <- F1
    cm_result[["MCC"]] <- MCC


    return(cm_result)

  }



  #' MLP_m
  #' @title Predict Interactions using Multi-layer Perception (MLP)
  #' @param data A matrix containing concatenated co-elution profiles,
  #' generated from \code{\link{getPPI}}.
  #' @param train_d A matrix of training data containing numerical features,
  #' generated from \code{\link{build_trainingData}}.
  #' @param train_l A vector of binary categorical label (0-1),
  #' generated from \code{\link{build_trainingData}}.
  #' @param nlayers Number of hidden layers. Defaults to 2.
  #' @param powerto1 Integer, the number of neurons in the first hidden layer
  #' as defined by two to the power of this value.Defaults to 6.
  #' @param powerto2 Integer, the number of neurons in the subsequent hidden
  #' layer as defined by two to the power of this value. Defaults to 7.
  #' @param drate Numeric, the dropout rates range. Defaults to 0.1.
  #' @param optimizer Name of the optimizee. For most models,
  #' this defaults to "rmsprop".
  #' @param b_size Number of samples per gradient update. Defaults to 128.
  #' @param epochs Number of epochs to train the model. Defaults to 50.
  #' @param cv_fold Number of partitions for cross-validation; defaults to 5.
  #' @param plots Logical value, indicating whether to plot the performance of
  #' the learning algorithm using k-fold cross-validation; defaults to FALSE.
  #' \itemize{ \item{pr_plot} - Precision-recall plot
  #' \item{roc_plot} - ROC plot
  #' \item{radar_plot} - Radar plot showing
  #' accuracy, F1-score , positive predictive value (PPV), sensitivity (SE)
  #' and MCC.}
  #' @param filename character string, indicating the output filename as an pdf
  #' object. Defaults to plots.pdf.
  #' @return Predicted interactions with predicted scores.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom dplyr rename
  #' @importFrom pROC roc
  #' @importFrom PRROC pr.curve
  #' @importFrom pROC ggroc
  #' @importFrom dplyr inner_join
  #' @importFrom stats predict
  #' @importFrom ggplot2 aes
  #' @importFrom ggplot2 ggplot
  #' @importFrom ggplot2 annotate
  #' @importFrom ggplot2 coord_equal
  #' @importFrom ggplot2 element_blank
  #' @importFrom ggplot2 element_line
  #' @importFrom ggplot2 element_text
  #' @importFrom ggplot2 geom_abline
  #' @importFrom ggplot2 scale_y_continuous
  #' @importFrom ggplot2 theme
  #' @importFrom ggplot2 theme_bw
  #' @importFrom ggplot2 unit
  #' @importFrom ggplot2 xlab
  #' @importFrom ggplot2 geom_line
  #' @importFrom ggplot2 ylab
  #' @importFrom grDevices dev.off
  #' @importFrom grDevices pdf
  #' @importFrom grDevices recordPlot
  #' @description This function uses the feedforward feedforward deep neural
  #' networks (DNN) (aka multi-layer perceptron)  to
  #' predict interactions from co-elution data.
  #' @export
  #' @examples
  #' #load the co-elution data
  #' data("HelaCE")
  #' #load the reference data for training
  #' data("refcpx")
  #' # concatenate the profile
  #' m_combined <- getPPI(HelaCE, similarity_calculation = TRUE)
  #' # build training data
  #' t_data <- build_trainingData(HelaCE, refcpx)
  #' #predict
  #' pred_int <-
  #' MLP_m(m_combined,
  #' t_data$train_d,
  #' t_data$train_l,
  #' cv_fold = 2)


  MLP_m <-
    function(data,
             train_d,
             train_l,
             nlayers = 2,
             powerto1 = 6,
             powerto2 = 5,
             drate = 0.1,
             optimizer = "rmsprop",
             epochs = 50,
             b_size = 128, cv_fold = 5, plots = FALSE,
             filename="plots.pdf"){

      if(nrow(train_d) != length(train_l)) {
        stop("Size of the class lables must be equal to the number of
             interactions in the training set")

      }

      if(!is.matrix(train_d)) train_d <- as.matrix(train_d)

      if(is.character(train_d) == TRUE){
        stop("Training matrix must include numerical variables")
      }

      if(is.character(train_l) == TRUE){
        stop("Class labels must be a numerical vector")
      }

     if(!is.vector(train_l)) train_l <- as.vector(train_l)


     if(is.character(data) == TRUE){
       stop("Data matrix must include numerical variables")
     }

     if(is.null(row.names(data))){
       stop("Please specify the row.names of co-eluting data")
     }

     . <- NULL
     V1 <- NULL
     V2 <- NULL

    pdata <- data
    PPI <- row.names(pdata)


    # define the  network
    network <- keras_model_sequential() %>%
      layer_dense(units = 2^powerto1, activation = "relu",
                  input_shape = c(ncol(train_d)),
                  kernel_regularizer = regularizer_l2(0.001)) %>%
      layer_batch_normalization() %>%
      layer_dropout(rate = drate)

    if(nlayers > 1) {
      map(2:nlayers, ~ network %>%
            layer_dense(units = 2^powerto2, activation = "relu",
                        kernel_regularizer = regularizer_l2(0.001)) %>%
            layer_batch_normalization() %>%
            layer_dropout(rate = drate)
      )

    }

    network %>%
      layer_dense(units = 1, activation = "sigmoid")

    network %>% compile(
      optimizer = optimizer,
      loss = "binary_crossentropy",
      metrics = c('accuracy'))


    # K-fold perfomance on training set and all data
    ptr <- pdata
    prob_u = matrix(nrow = nrow(ptr), ncol = 1)
    index_u = rep(1:cv_fold, nrow(ptr))
    ind_u = index_u[1:nrow(ptr)]



    xtr <- train_d
    ytr <- train_l
    prob = matrix(nrow = length(ytr), ncol = 1)
    index = rep(1:cv_fold, nrow(xtr))
    ind = index[1:nrow(xtr)]


    for (k in 1:cv_fold) {
      cat(".")

      # training data
      xcal <- xtr[ind != k, ]
      ycal <- ytr[ind != k]
      xtest <- xtr[ind == k, ]

      # all data
      xtest_u <- ptr[ind_u == k, ]

      # train the network
      history <- network %>% fit(
        xcal,
        ycal,
        verbose = 0,
        epochs = epochs,
        batch_size = b_size,
        callbacks =
          list(callback_early_stopping(monitor = "loss",patience = 5,
                                       restore_best_weights = TRUE)))

      prob[ind == k, ] =  network %>% predict(xtest)
      prob_u[ind_u == k, ] =  network %>% predict(xtest_u)

    }

    prob_u <-
      as.data.frame(prob_u) %>%
      cbind(PPI, .)
    colnames(prob_u) <- c("PPI", "Prob")



    # add plot here
    if (plots) {

      # Performance evaluation
      prob <- as.data.frame(prob)
      colnames(prob) <-c("Prob")
      categ_lab <- ifelse(prob$Prob >= 0.5, 1, 0)


      cm_out <- .cm_f(categ_lab , ytr)
      # roc analysis
      roc_c <-
        roc(ytr, prob$Prob)
      # pr analysis
      pr_c <-
        pr.curve(
          scores.class0 = prob$Prob[ytr == 1],
          scores.class1 =  prob$Prob[ytr == 0],
          curve = TRUE)


      pdf(filename)

      # Generate plot for cm result
      df <-
        data.frame(rbind(rep(1,6),
                         rep(0,6), cbind(cm_out$ACC,cm_out$SE,cm_out$SP,
                                         cm_out$PPV,cm_out$F1,cm_out$MCC)))
      colnames(df) <- c("ACC","SE","SP","PPV","F1", "MCC")

      raderplot <-
        fmsb::radarchart(df,cglty = 2, pfcol = c("#99999980"),
                   cglcol = "blue",pcol = 2,plwd = 2, plty = 1)

      p <- recordPlot()



      # roc curve
      roc_plot <-
        ggroc(roc_c, legacy.axes = TRUE) +
        geom_abline(
          slope = 1, intercept = 0,
          linetype = "dashed", alpha = 0.7, color = "darkgreen"
        ) + coord_equal() + theme_bw() +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        theme(axis.ticks.length = unit(.5, "cm")) +
        theme(text = element_text(size = 14, color = "black")) +
        theme(axis.text = element_text(size = 12, color = "black")) +
        xlab("False Positive Rate (1-Specificity)") +
        ylab("True Positive Rate (Sensitivity)") +
        annotate("text", x=0.25, y=0.8,
                 label= paste("AUC:",round(roc_c[["auc"]],2)))


      # PR_curve
      PR_Object <- as.data.frame(pr_c$curve)

      pr_plot <-
        ggplot(PR_Object, aes(x = V1, y = V2)) +
        geom_line() +
        theme_bw() + scale_y_continuous(limits = c(0, 1)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
        theme(axis.ticks.length = unit(.5, "cm")) +
        theme(text = element_text(size = 14, color = "black")) +
        theme(axis.text = element_text(size = 12, color = "black")) +
        xlab("Recall") + ylab("Percision") +
        annotate("text", x=0.25, y=0.25,
                 label= paste("AUC:",round(pr_c[["auc.integral"]],2)))

      # Make plots.
      plot_list <- list()
      plot_list$roc <- roc_plot
      plot_list$pr <- pr_plot
      plot_list$rader <- p


      for (i in seq_along(plot_list)) {
        print(plot_list[[i]])
      }
      dev.off()

    }


    return(prob_u)

    }












